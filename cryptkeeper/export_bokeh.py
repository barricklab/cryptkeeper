# Visualize the data on a bokah plot
from bokeh.models import ColumnDataSource, Grid, LinearAxis, Plot, Rect, TextInput, CustomJS
from bokeh.io import curdoc, show
import numpy as np
from collections import namedtuple
import pandas as pd
from bokeh.transform import linear_cmap, log_cmap
from bokeh.palettes import Viridis256
from bokeh.plotting import output_file, save
from bokeh.layouts import column, row

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


def export_bokeh(cryptresult, filename=None):
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

    wigits = []

    fig.extra_y_ranges = {"y_range2": Range1d(start=-5000, end=15000),
                          "y_range3": Range1d(start=0, end=0)}
    fig.add_layout(LinearAxis(y_range_name="y_range3"), 'left')

    # Set size of the figure to the width of the browser
    fig.width = 1000

    # Draw a line for the sequence
    fig.line([0, len(cryptresult.sequence)], [0, 0], line_width=2, color="black", y_range_name="y_range2")

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
            annotation_base_y = -1250.*genbank_annotation.nest_level-3000
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

        genbank_glyphs = fig.patches('x', 'y', color='color', source=genbank_dictionary, alpha=0.5, line_color='black', line_width=1, y_range_name="y_range2")
        genbank_glyphs_hover = HoverTool(renderers=[genbank_glyphs], tooltips=[("Name", "@name")])
        fig.add_tools(genbank_glyphs_hover)

    # Add the expressed CDSs
    expressed_CDSs = cryptresult.translation_sites
    if expressed_CDSs:
        # Sort the expressed_CDS by the difference between the start and end
        def sort_algorythm(x):
            candidate_1 = x.end - x.start
            candidate_2 = x.start - x.end
            return max(candidate_1, candidate_2)


        expressed_CDSs = sorted(expressed_CDSs, key=sort_algorythm, reverse=True)

        boxes, highest_expression = plot_boxes(expressed_CDSs)
        source = ColumnDataSource(boxes)

        cmap = linear_cmap(field_name='burden', palette=Viridis256, low=0, high=highest_expression)
        rectangles = fig.rect(x="x", y="y", width="w", height="h", source=source, color=cmap, line_color="black", line_width=1, y_range_name="y_range3")
        rectangles_hover = HoverTool(renderers=[rectangles], tooltips=[("Burden", "@burden")])
        fig.add_tools(rectangles_hover)

        color_bar = rectangles.construct_color_bar(padding=0,
                                        ticker=fig.xaxis.ticker,
                                        formatter=fig.xaxis.formatter)
        fig.add_layout(color_bar, 'right')
        highest_y = max([(x[0]+x[1]) for x in zip(boxes['y'], boxes['h'])])

        max_y = TextInput(title="Max y", value=str(highest_y))
        max_y_js = CustomJS(args=dict(plotRange1=fig.y_range, plotRange2=fig.extra_y_ranges["y_range3"]), code="""
    var newEnd = Number(cb_obj.value)
    plotRange2.start = newEnd * -1/3     
    plotRange2.end = newEnd                        
""")
        wigits.append((max_y, max_y_js))

        max_burden = TextInput(title="Max burden", value=str(highest_expression))
        max_burden_js = CustomJS(args=dict(color_bar=color_bar), code="""
    var newMax = Number(cb_obj.value)
    color_bar.color_mapper.high = newMax
""")
        wigits.append((max_burden, max_burden_js))

    # add promoters
    promoters = cryptresult.promoters
    if promoters:
        promoters = sorted(promoters, key=lambda x: x.score, reverse=True)
        promoter_dict = {'x': [],
                     'y': [],
                     'position': [],
                     'score': []}
        
        for promoter in promoters[:10]:
            arrow_shape = ((-10, 0), (10, 0), (10, 400), (50, 400), (50, 350), (100, 450), (50, 550), (50, 500), (-10, 500), (-10, 0))
            if promoter.strand == "-":
                arrow_shape = [(-x[0], -x[1]) for x in arrow_shape]
            xs = [promoter.TSSpos + x[0] for x in arrow_shape]
            ys = [x[1]-1000 for x in arrow_shape]
            promoter_dict['x'].append(xs)
            promoter_dict['y'].append(ys)
            promoter_dict['position'].append(promoter.TSSpos)
            promoter_dict['score'].append(promoter.score)


        promoter_glyphs = fig.patches('x', 'y', color='green', source=promoter_dict, alpha=0.5, line_color='black', line_width=1, y_range_name="y_range2")
        promoter_glyphs_hover = HoverTool(renderers=[promoter_glyphs], tooltips=[("Position", "@position"), ("Score", "@score")])
        fig.add_tools(promoter_glyphs_hover)

    # Add the terminators
    rdpt = cryptresult.rho_dep_terminators
    ridpt = cryptresult.rho_ind_terminators

    if rdpt:
        rdpt = sorted(rdpt, key=lambda x: x.scores[0], reverse=True)
        t_shape = ((-10, 0), (10, 0), (10, 400), (50, 400), (50, 500), (-50, 500), (-50, 400), (-10, 400), (-10, 0))
        terminator_dict = {'x': [],
                           'y': [],
                           'position': [],
                           'score': []}
        # strand, start_rut, end_rut,  score
        for terminator in rdpt[0:5]:
            if terminator.strand == "-":
                t_shape = [(-x[0], -x[1]) for x in t_shape]
            xs = [terminator.start_rut + x[0] for x in t_shape]
            ys = [x[1]-1000 for x in t_shape]
            terminator_dict['x'].append(xs)
            terminator_dict['y'].append(ys)
            terminator_dict['position'].append(f'{terminator.start_rut}-{terminator.end_rut}')
            terminator_dict['score'].append(terminator.scores[0])

        terminator_glyphs = fig.patches('x', 'y', color='red', source=terminator_dict, alpha=0.5, line_color='black', line_width=1, y_range_name="y_range2")
        terminator_glyphs_hover = HoverTool(renderers=[terminator_glyphs], tooltips=[("Position", "@position"), ("Score", "@score")])
        fig.add_tools(terminator_glyphs_hover)

    if ridpt:
        ridpt = sorted(ridpt, key=lambda x: x.conf, reverse=True)
        t_shape = ((-10, 0), (10, 0), (10, 400), (50, 400), (50, 500), (-50, 500), (-50, 400), (-10, 400), (-10, 0))
        terminator_dict = {'x': [],
                            'y': [],
                            'position': [],
                            'score': []}
        for terminator in ridpt[0:5]:
            # strand, conf, start, end
            if terminator.strand == "-":
                t_shape = [(-x[0], -x[1]) for x in t_shape]
            xs = [terminator.start + x[0] for x in t_shape]
            ys = [x[1]-1000 for x in t_shape]
            terminator_dict['x'].append(xs)
            terminator_dict['y'].append(ys)
            terminator_dict['position'].append(f'{terminator.start}-{terminator.end}')
            terminator_dict['score'].append(terminator.conf)

        terminator_glyphs = fig.patches('x', 'y', color='blue', source=terminator_dict, alpha=0.5, line_color='black', line_width=1, y_range_name="y_range2")
        terminator_glyphs_hover = HoverTool(renderers=[terminator_glyphs], tooltips=[("Position", "@position"), ("Score", "@score")])
        fig.add_tools(terminator_glyphs_hover)

    # Extra line below promoters and terminators
    fig.line([0, len(cryptresult.sequence)], [-2000, -2000], line_width=2, color="black", y_range_name="y_range2")

    # set y axis min and max of y_range3
    fig.extra_y_ranges["y_range3"].start = -highest_y/3
    fig.extra_y_ranges["y_range3"].end = highest_y

    widgets_to_add = []
    for wigit, js in wigits:
        wigit.js_on_change('value', js)
        widgets_to_add.append(wigit)

    # Add the plot to the document and add a title
    layout = column(fig, *widgets_to_add)
    layout.sizing_mode = 'scale_width'


    # show the results
    show(layout)

    if filename:
        output_file(filename=filename, title="Static HTML file")
        save(fig)


    return fig


# EOF
