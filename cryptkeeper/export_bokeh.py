# Visualize the data on a bokah plot
import math
from bokeh.models import ColumnDataSource, Grid, LinearAxis, Plot, Rect
from bokeh.io import curdoc, show
import numpy as np
from collections import namedtuple
import pandas as pd
from bokeh.transform import linear_cmap, log_cmap
from bokeh.palettes import Viridis256
from bokeh.plotting import output_file, save

def plot_boxes(features_list):

    placed_ORFs = pd.DataFrame(columns=["start_x", "stop_x", "start_y", "stop_y", "burden"])
    highest_expression = 0
    highest_burden = 0

    for feature in features_list:
        start_x_adding = int(feature.start)
        stop_x_adding = int(feature.end)
        burden = float(feature.burden)
        if burden > highest_burden:
            highest_burden = burden

        
        if start_x_adding > stop_x_adding:  # Then the ORF is reverse stranded
            start_x_adding, stop_x_adding = stop_x_adding, start_x_adding
            # Now the variables are switched
        
        expression = float(10 ** abs(feature.expression))
        if expression > highest_expression:
            highest_expression = expression


        to_add = 0

        failed = True
        while failed:
            start_y_adding = to_add
            stop_y_adding = to_add + expression  # Subtracts if reverse strand
            about_to_break = False

            for j in range(placed_ORFs.shape[0]):
                start_x_checking = int(placed_ORFs.loc[j, "start_x"])
                stop_x_checking = int(placed_ORFs.loc[j, "stop_x"])
                start_y_checking = float(placed_ORFs.loc[j, "start_y"])
                stop_y_checking = float(placed_ORFs.loc[j, "stop_y"])

                if start_x_adding >= stop_x_checking:
                    continue
                if stop_x_adding <= start_x_checking:
                    continue
                
                if np.sign(expression) != np.sign(stop_y_checking):
                    continue
                
                if np.sign(expression) == 1:
                    if start_y_adding >= stop_y_checking:
                        continue
                    if stop_y_adding <= start_y_checking:
                        continue
                    to_add = stop_y_checking
                    about_to_break = True
                    break
                
                if np.sign(expression) == -1:
                    if start_y_adding <= stop_y_checking:
                        continue
                    if stop_y_adding >= start_y_checking:
                        continue
                    to_add = stop_y_checking
                    about_to_break = True
                    break
                
                about_to_break = True
                break
            
            if not about_to_break:
                failed = False
        
        placed_ORFs.loc[placed_ORFs.shape[0]] = [start_x_adding, stop_x_adding, start_y_adding, stop_y_adding, burden]


    # Convert to somethong bokah understands
    bokah_orfs = {
        "x": [],
        "y": [],
        "w": [],
        "h": [],
        "burden": [],
    }

    for i in range(placed_ORFs.shape[0]):
        bokah_orfs["x"].append((placed_ORFs.loc[i, "start_x"] + placed_ORFs.loc[i, "stop_x"]) / 2)
        bokah_orfs["y"].append((placed_ORFs.loc[i, "start_y"] + placed_ORFs.loc[i, "stop_y"]) / 2)
        bokah_orfs["w"].append(placed_ORFs.loc[i, "stop_x"] - placed_ORFs.loc[i, "start_x"])
        bokah_orfs["h"].append(abs(placed_ORFs.loc[i, "stop_y"] - placed_ORFs.loc[i, "start_y"]))
        bokah_orfs["burden"].append(placed_ORFs.loc[i, "burden"])




    return bokah_orfs, highest_burden


def export_bokah(cryptresult, filename=None):
    # Set up plot
    from bokeh.plotting import figure, show
    from bokeh.models import ColumnDataSource, Range1d, LabelSet, Label
    from bokeh.models import Arrow, OpenHead, NormalHead, VeeHead
    from bokeh.models import BoxAnnotation
    from bokeh.models import HoverTool
    from bokeh.models import Legend
    from bokeh.models import Span
    from bokeh.models import Arrow, OpenHead, NormalHead, VeeHead
    from bokeh.models import BoxAnnotation
    from bokeh.models import HoverTool
    from bokeh.models import Legend
    from bokeh.models import Span
    from bokeh.models import PolyAnnotation

    # Set up the figure
    fig = figure()
    fig.xaxis.axis_label = "Position"
    fig.yaxis.axis_label = "Strand"
    fig.yaxis.visible = False
    fig.ygrid.visible = False
    fig.xgrid.visible = False
    fig.toolbar.logo = None

    # Set size of the figure to the width of the browser
    fig.width = 1000

    # Draw a line for the sequence
    fig.line([0, len(cryptresult.sequence)], [0, 0], line_width=2, color="black")

    # Add predefined features from the genbank file
    if cryptresult.annotations:
        from Bio import SeqIO
        from Bio.SeqFeature import FeatureLocation
        from Bio.SeqFeature import CompoundLocation
        from Bio.SeqFeature import ExactPosition
        from Bio.SeqFeature import BeforePosition
        from Bio.SeqFeature import AfterPosition
        from Bio.SeqFeature import UnknownPosition
        from Bio.SeqFeature import FeatureLocation
        from Bio.SeqFeature import CompoundLocation
        from Bio.SeqFeature import ExactPosition
        from Bio.SeqFeature import BeforePosition
        from Bio.SeqFeature import AfterPosition
        from Bio.SeqFeature import UnknownPosition

        genbank_dictionary = {'x': [],
                              'y': [],
                              'color': [],
                              'name': [],
                              }

        
        for genbank_annotation in cryptresult.annotations:

            arrow_depth = 100
            if genbank_annotation.end-genbank_annotation.start < arrow_depth:
                    arrow_depth = genbank_annotation.end-genbank_annotation.start
            annotation_base_y = -2000.*genbank_annotation.nest_level-1000
            if genbank_annotation.strand == 1:
                xs=[genbank_annotation.start, genbank_annotation.start, genbank_annotation.end-arrow_depth, genbank_annotation.end, genbank_annotation.end-arrow_depth]
                ys=[annotation_base_y+500, annotation_base_y-500, annotation_base_y-500, annotation_base_y, annotation_base_y+500]
            elif genbank_annotation.strand == -1:
                xs=[genbank_annotation.end, genbank_annotation.end, genbank_annotation.start+arrow_depth, genbank_annotation.start, genbank_annotation.start+arrow_depth]
                ys=[annotation_base_y-500, annotation_base_y+500, annotation_base_y+500, annotation_base_y, annotation_base_y-500]
            else:
                xs=[genbank_annotation.start, genbank_annotation.start, genbank_annotation.end, genbank_annotation.end]
                ys=[annotation_base_y-500, annotation_base_y+500, annotation_base_y+500, annotation_base_y-500]
            name = genbank_annotation.name
            color = genbank_annotation.color

            genbank_dictionary['x'].append(xs)
            genbank_dictionary['y'].append(ys)
            genbank_dictionary['color'].append(color)
            genbank_dictionary['name'].append(name)

        genbank_glyphs = fig.patches('x', 'y', color='color', source=genbank_dictionary, alpha=0.5, line_color='black', line_width=1)
        genbank_glyphs_hover = HoverTool(renderers=[genbank_glyphs], tooltips=[("Name", "@name")])
        fig.add_tools(genbank_glyphs_hover)



        source = ColumnDataSource(genbank_dictionary)

        # Add mouseover text from feature name
        #hover = HoverTool(tooltips=[("Name", "@name")])
        #fig.add_tools(hover)

        # Draw the features
        #fig.rect(x="x", y="y", width="w", height="h", source=source, color="color", line_color="black", line_width=1)
        

    # Add the expressed CDSs
    expressed_CDSs = cryptresult.translation_sites

    # Sort the expressed_CDS by the difference between the start and end
    def sort_algorythm(x):
        candidate_1 = x.end - x.start
        candidate_2 = x.start - x.end
        return max(candidate_1, candidate_2)


    expressed_CDSs = sorted(expressed_CDSs, key=sort_algorythm, reverse=True)

    boxes, highest_expression = plot_boxes(expressed_CDSs)
    source = ColumnDataSource(boxes)

    cmap = linear_cmap(field_name='burden', palette=Viridis256, low=0, high=highest_expression)
    rectangles = fig.rect(x="x", y="y", width="w", height="h", source=source, color=cmap, line_color="black", line_width=1)
    rectangles_hover = HoverTool(renderers=[rectangles], tooltips=[("Burden", "@burden")])
    fig.add_tools(rectangles_hover)

    color_bar = rectangles.construct_color_bar(padding=0,
                                      ticker=fig.xaxis.ticker,
                                      formatter=fig.xaxis.formatter)
    fig.add_layout(color_bar, 'right')


    # show the results
    show(fig)

    if filename:
        output_file(filename=filename, title="Static HTML file")
        save(fig)


    return fig


# EOF
