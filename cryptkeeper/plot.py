# Visualize the data on a bokah plot

import os
import copy
import shutil
import math
import bokeh
import numpy as np
import pandas as pd
from bokeh.transform import linear_cmap
from bokeh.palettes import viridis
from bokeh.layouts import column, row
from bokeh.events import DocumentReady
from bokeh.io import curdoc, export_svg
from bokeh.plotting import figure
from bokeh.models import (
    ColumnDataSource,
    Range1d,
    HoverTool,
    LinearAxis,
    TextInput,
    CustomJS,
    Div,
    ResetTool,
    WheelZoomTool,
    PanTool,
    BoxZoomTool,
    SaveTool,
    NumberFormatter,
)
from bokeh.embed import components
from bokeh.models.widgets import DataTable, TableColumn
from bokeh.core.validation import silence
from bokeh.core.validation.warnings import MISSING_RENDERERS

import logging

logging.getLogger("bokeh.io.export").setLevel(logging.CRITICAL)

ANNOTATION_SPACE = 25  # In percent of the graph
GLYPH_SCALING = 1.2
GLYPHS_PER_KB = 3

FONT = "Arial"
FONTSIZE = "12pt"


def plot_boxes(features_list):
    """
    Generate the plot of boxes for a list of features.
    Parameters:
        features_list (list): A list of features, typically the expressed ORFs.
    Returns:
        tuple: A tuple containing the dictionary of bokah ORFs with their coordinates and attributes, and the highest burden value.
    """
    placed_ORFs = pd.DataFrame(
        columns=[
            "start_x",
            "stop_x",
            "start_y",
            "stop_y",
            "burden",
            "strand",
            "tir",
        ]
    )
    highest_tir = 0
    highest_burden = 0

    for feature in features_list:
        start_x = int(feature.start)
        stop_x = int(feature.end)
        burden = float(feature.burden)
        strand = feature.strand
        tir = feature.tir

        if burden > highest_burden:
            highest_burden = burden

        if start_x > stop_x:  # Reverse stranded ORF
            start_x, stop_x = stop_x, start_x

        # tir = float(10 ** abs(tir))  # Calculate tir
        if tir > highest_tir:
            highest_tir = tir

        to_add = 0
        failed = True

        start_y = 0
        stop_y = 0 + tir

        while failed:
            start_y = to_add
            stop_y = to_add + tir

            about_to_break = False

            for j in range(placed_ORFs.shape[0]):
                start_x_checking = int(placed_ORFs.loc[j, "start_x"])
                stop_x_checking = int(placed_ORFs.loc[j, "stop_x"])
                start_y_checking = float(placed_ORFs.loc[j, "start_y"])
                stop_y_checking = float(placed_ORFs.loc[j, "stop_y"])

                if start_x >= stop_x_checking or stop_x <= start_x_checking:
                    continue

                if np.sign(tir) != np.sign(stop_y_checking):
                    continue

                if np.sign(tir) == 1:
                    if start_y >= stop_y_checking or stop_y <= start_y_checking:
                        continue
                    to_add = stop_y_checking
                    about_to_break = True
                    break

                if np.sign(tir) == -1:
                    if start_y <= stop_y_checking or stop_y >= start_y_checking:
                        continue
                    to_add = stop_y_checking
                    about_to_break = True
                    break

                about_to_break = True
                break

            if not about_to_break:
                failed = False

        placed_ORFs.loc[placed_ORFs.shape[0]] = [
            start_x,
            stop_x,
            start_y,
            stop_y,
            burden,
            strand,
            tir,
        ]

    bokah_orfs = {
        "x": [],
        "y": [],
        "w": [],
        "h": [],
        "rbs_strength": [],
        "burden": [],
        "position": [],
        "strand": [],
    }

    for i in range(placed_ORFs.shape[0]):
        bokah_orfs["x"].append(
            (placed_ORFs.loc[i, "start_x"] + placed_ORFs.loc[i, "stop_x"]) / 2
        )
        bokah_orfs["y"].append(
            (placed_ORFs.loc[i, "start_y"] + placed_ORFs.loc[i, "stop_y"]) / 2
        )
        bokah_orfs["w"].append(
            placed_ORFs.loc[i, "stop_x"] - placed_ORFs.loc[i, "start_x"]
        )
        bokah_orfs["h"].append(
            abs(placed_ORFs.loc[i, "stop_y"] - placed_ORFs.loc[i, "start_y"])
        )
        bokah_orfs["rbs_strength"].append(placed_ORFs.loc[i, "tir"])
        bokah_orfs["burden"].append(placed_ORFs.loc[i, "burden"])
        bokah_orfs["position"].append(
            f'{placed_ORFs.loc[i, "start_x"]}-{placed_ORFs.loc[i, "stop_x"]}'
        )
        bokah_orfs["strand"].append(f'{placed_ORFs.loc[i, "strand"]}')

    return bokah_orfs, highest_burden


def make_plot(cryptresult, tick_frequency=1000, filename=None, show_small=False):
    view_format = "mirrored"
    # Set up plot

    annotation_scale_range = len(cryptresult.sequence)
    number_of_glyphs = math.floor(len(cryptresult.sequence) / 1000 * GLYPHS_PER_KB)

    # Set up the figure
    fig = figure(
        width=1500,
        height=750,
        tools=[ResetTool(), WheelZoomTool(), PanTool(), BoxZoomTool(), SaveTool()],
    )

    fig.xaxis.axis_label = "Position"
    fig.yaxis.axis_label = "Predicted Translation"
    fig.yaxis.visible = False
    fig.ygrid.visible = False
    fig.xgrid.visible = False
    fig.toolbar.logo = None

    fig.lod_threshold = None

    wigits = []
    tables = {}

    fig.extra_y_ranges = {
        "y_range2": Range1d(start=-1, end=1),
        "y_range3": Range1d(start=-1, end=1),
        "y_range4": Range1d(start=-1, end=1),
    }
    shownaxis = LinearAxis(y_range_name="y_range3", axis_label="Predicted Translation")
    fig.add_layout(shownaxis, "left")

    fig.extra_x_ranges = {
        "x_range3": Range1d(start=0, end=len(cryptresult.sequence)),
    }

    fig.xaxis.axis_label_text_font_size = FONTSIZE
    fig.yaxis.axis_label_text_font_size = FONTSIZE
    shownaxis.axis_label_text_font_size = FONTSIZE
    fig.xaxis.axis_label_text_font = FONT
    fig.yaxis.axis_label_text_font = FONT
    shownaxis.axis_label_text_font = FONT
    fig.xaxis.axis_label_text_font_style = "bold"
    fig.yaxis.axis_label_text_font_style = "bold"
    shownaxis.axis_label_text_font_style = "bold"
    fig.xaxis.major_label_text_font_size = FONTSIZE
    fig.yaxis.major_label_text_font_size = FONTSIZE
    shownaxis.major_label_text_font_size = FONTSIZE
    fig.xaxis.major_label_text_font = FONT
    fig.yaxis.major_label_text_font = FONT
    shownaxis.major_label_text_font = FONT

    # set width to length of sequence
    fig.x_range = Range1d(start=0, end=len(cryptresult.sequence))

    # Set up the annotation tracks
    forward_exists = False
    reverse_exists = False
    for feature in (
        cryptresult.rho_dep_terminators
        + cryptresult.int_terminators
        + cryptresult.promoters
    ):
        if feature.strand == "+" and not forward_exists:
            forward_exists = True
        elif feature.strand == "-" and not reverse_exists:
            reverse_exists = True
        elif forward_exists and reverse_exists:
            break

    # Draw a line for the top of the annotations
    fig.line(
        [0, annotation_scale_range],
        [0, 0],
        line_width=2,
        color="black",
        y_range_name="y_range2",
    )

    reverse_height = None
    forward_height = None
    if view_format == "stacked" and reverse_exists:
        forward_exists = True  # The reverse track goes in place of the forward track
        forward_height = -1000
        reverse_height = -1000
        fig.line(
            [0, annotation_scale_range],
            [-2000, -2000],
            line_width=2,
            color="black",
            y_range_name="y_range2",
        )
        annotation_depth = -3000
    elif forward_exists:
        # we add the track for forward strand annotations. We will add reverse height later
        forward_height = -1000
        fig.line(
            [0, annotation_scale_range],
            [-1500, -1500],
            line_width=2,
            color="black",
            y_range_name="y_range2",
        )
        annotation_depth = -2500
    elif view_format == "mirrored" and reverse_exists:
        annotation_depth = -500  # We will add the reverse track later after we know how much space annotations take up
    else:  # if neither exist
        annotation_depth = -500
    # Now we add the annotations so we know how much space it takes

    lowest_annotation_y = copy.copy(annotation_depth) - 500
    if cryptresult.annotations:
        # annotation_depth -= 500
        genbank_dictionary = {
            "x": [],
            "y": [],
            "color": [],
            "name": [],
            "position": [],
            "strand": [],
        }

        for genbank_annotation in cryptresult.annotations:
            arrow_depth = 100
            # @TODO: SCALE ANNOTATIONS

            scaled_start = (
                genbank_annotation.start
                * len(cryptresult.sequence)
                / annotation_scale_range
            )
            scaled_end = (
                genbank_annotation.end
                * len(cryptresult.sequence)
                / annotation_scale_range
            )

            if genbank_annotation.end - genbank_annotation.start < arrow_depth:
                arrow_depth = genbank_annotation.end - genbank_annotation.start
            annotation_base_y = (
                -1250.0 * genbank_annotation.nest_level + annotation_depth
            )

            if genbank_annotation.strand == 1:
                xs = [
                    scaled_start,
                    scaled_start,
                    scaled_end - arrow_depth,
                    scaled_end,
                    scaled_end - arrow_depth,
                ]
                ys = [
                    annotation_base_y + 500,
                    annotation_base_y - 500,
                    annotation_base_y - 500,
                    annotation_base_y,
                    annotation_base_y + 500,
                ]
            elif genbank_annotation.strand == -1:
                xs = [
                    scaled_end,
                    scaled_end,
                    scaled_start + arrow_depth,
                    scaled_start,
                    scaled_start + arrow_depth,
                ]
                ys = [
                    annotation_base_y - 500,
                    annotation_base_y + 500,
                    annotation_base_y + 500,
                    annotation_base_y,
                    annotation_base_y - 500,
                ]
            else:
                xs = [scaled_start, scaled_start, scaled_end, scaled_end]
                ys = [
                    annotation_base_y - 500,
                    annotation_base_y + 500,
                    annotation_base_y + 500,
                    annotation_base_y - 500,
                ]
            name = genbank_annotation.name

            # Update lowest annotation y
            if annotation_base_y - 500 < lowest_annotation_y:
                lowest_annotation_y = annotation_base_y - 500

            # Define the color
            color = genbank_annotation.color

            # Add the annotation to the dictionary

            genbank_dictionary["x"].append(xs)
            genbank_dictionary["y"].append(ys)
            genbank_dictionary["color"].append(color)
            genbank_dictionary["name"].append(name)
            genbank_dictionary["position"].append(
                f"{genbank_annotation.start}-{genbank_annotation.end}"
            )
            genbank_dictionary["strand"].append(genbank_annotation.strand)

        genbank_glyphs = fig.patches(
            "x",
            "y",
            color="color",
            source=genbank_dictionary,
            alpha=0.5,
            line_color="black",
            line_width=1,
            y_range_name="y_range2",
        )
        genbank_glyphs_hover = HoverTool(
            renderers=[genbank_glyphs],
            tooltips=[("Type", "Annotation"), ("Name", "@name")],
            visible=False,
        )
        fig.add_tools(genbank_glyphs_hover)

        # Draw a line below the annotations

    lowest_annotation_y -= 500
    fig.line(
        [0, annotation_scale_range],
        [lowest_annotation_y, lowest_annotation_y],
        line_width=2,
        color="black",
        y_range_name="y_range2",
    )
    annotation_depth = lowest_annotation_y

    if view_format == "mirrored" and reverse_exists:
        # We add the track for reverse strand annotations
        annotation_depth -= 500
        reverse_height = lowest_annotation_y - 500
        fig.line(
            [0, annotation_scale_range],
            [reverse_height - 1000, reverse_height - 1000],
            line_width=2,
            color="black",
            y_range_name="y_range2",
        )
        annotation_depth = reverse_height - 1000

    # add promoters
    promoters = cryptresult.promoters
    if promoters:
        promoters = sorted(promoters, key=lambda x: x.score, reverse=True)
        promoter_dict = {"x": [], "y": [], "position": [], "score": [], "strand": []}

        for promoter in promoters:
            arrow_shape = (
                (0, 0),
                (20, 0),
                (20, 350),
                (40, 350),
                (40, 300),
                (90, 400),
                (40, 500),
                (40, 450),
                (0, 450),
                (0, 0),
            )
            scaled_arrow = [
                (
                    x[0] * len(cryptresult.sequence) / 10000 * GLYPH_SCALING,
                    x[1] * GLYPH_SCALING,
                )
                for x in arrow_shape
            ]
            arrow_shape = scaled_arrow
            scaled_position = promoter.TSSpos
            if promoter.strand == "-":
                arrow_shape = [(-x[0], -x[1]) for x in arrow_shape]
            xs = [scaled_position + x[0] for x in arrow_shape]
            if promoter.strand == "+":
                ys = [x[1] + forward_height for x in arrow_shape]
            else:
                ys = [x[1] + reverse_height for x in arrow_shape]
            promoter_dict["x"].append(xs)
            promoter_dict["y"].append(ys)
            promoter_dict["position"].append(promoter.TSSpos)
            promoter_dict["score"].append(promoter.score)
            promoter_dict["strand"].append(promoter.strand)

        promoter_glyphs = fig.patches(
            "x",
            "y",
            color="green",
            source=promoter_dict,
            alpha=0.5,
            line_color="black",
            line_width=1,
            y_range_name="y_range2",
        )
        promoter_glyphs_hover = HoverTool(
            renderers=[promoter_glyphs],
            tooltips=[
                ("Position", "@position"),
                ("Strand", "@strand"),
                ("Score", "@score"),
            ],
            visible=False,
        )
        fig.add_tools(promoter_glyphs_hover)

        # Add a widget to change the number of promoters shown
        promoter_number = TextInput(
            title="Number of promoters",
            value=str(number_of_glyphs),
        )
        promoter_javascript = """
                            var originalPromoData = structuredClone(source)
                            var newNumber = Number(promoter_number.value)
                            promoter_glyphs.data_source.data['x'] = originalPromoData['x'].slice(0, newNumber)
                            promoter_glyphs.data_source.data['y'] = originalPromoData['y'].slice(0, newNumber)
                            promoter_glyphs.data_source.data['position'] = originalPromoData['position'].slice(0, newNumber)
                            promoter_glyphs.data_source.data['score'] = originalPromoData['score'].slice(0, newNumber)
                            promoter_glyphs.data_source.data['strand'] = originalPromoData['strand'].slice(0, newNumber)
                            promoter_glyphs.data_source.change.emit()
                            """
        promoter_javascript = CustomJS(
            args=dict(
                promoter_glyphs=promoter_glyphs,
                promoter_number=promoter_number,
                source=promoter_dict,
            ),
            code=promoter_javascript,
        )

        curdoc().js_on_event(DocumentReady, promoter_javascript)
        fig.js_on_event("reset", promoter_javascript)
        wigits.append((promoter_number, promoter_javascript))

        # Add table
        name = "Promoters"
        promoter_table = generate_bokeh_table(promoters, name)
        tables[name] = promoter_table

    # Add the terminators
    rdpt = cryptresult.rho_dep_terminators
    ridpt = cryptresult.int_terminators

    if rdpt:
        rdpt = sorted(rdpt, key=lambda x: x.score, reverse=True)

        ritterminator_dict = {
            "x": [],
            "y": [],
            "position": [],
            "score": [],
            "strand": [],
        }
        # strand, start_rut, end_rut,  score
        for terminator in rdpt:
            t_shape = (
                (0, 0),
                (20, 0),
                (20, 350),
                (60, 350),
                (60, 450),
                (-40, 450),
                (-40, 350),
                (0, 350),
                (0, 0),
            )
            scaled_t = [
                (
                    x[0] * len(cryptresult.sequence) / 10000 * GLYPH_SCALING,
                    x[1] * GLYPH_SCALING,
                )
                for x in t_shape
            ]
            t_shape = scaled_t
            scaled_position = terminator.start_rut
            if terminator.strand == "-":
                t_shape = [(-x[0], -x[1]) for x in t_shape]
            xs = [scaled_position + x[0] for x in t_shape]
            if terminator.strand == "+":
                ys = [x[1] + forward_height for x in t_shape]
            else:
                ys = [x[1] + reverse_height for x in t_shape]
            ritterminator_dict["x"].append(xs)
            ritterminator_dict["y"].append(ys)
            ritterminator_dict["position"].append(
                f"{terminator.start_rut}-{terminator.end_rut}"
            )
            ritterminator_dict["score"].append(terminator.score)
            ritterminator_dict["strand"].append(terminator.strand)

        terminator_glyphs = fig.patches(
            "x",
            "y",
            color="red",
            source=ritterminator_dict,
            alpha=0.5,
            line_color="black",
            line_width=1,
            y_range_name="y_range2",
        )
        terminator_glyphs_hover = HoverTool(
            renderers=[terminator_glyphs],
            tooltips=[
                ("Type", "Rho-dependent Terminator"),
                ("Position", "@position"),
                ("Strand", "@strand"),
                ("Score", "@score"),
            ],
            visible=False,
        )
        fig.add_tools(terminator_glyphs_hover)

        # Add a widget to change the number of terminators shown
        rdpt_number = TextInput(
            title="Number of Rho-Dependant Terminators", value=str(number_of_glyphs)
        )
        rdpt_javascript = """
                        var originalRDPTData = structuredClone(source)
                        var newNumber = Number(rdpt_number.value)
                        rdpt_glyphs.data_source.data['x'] = originalRDPTData['x'].slice(0, newNumber)
                        rdpt_glyphs.data_source.data['y'] = originalRDPTData['y'].slice(0, newNumber)
                        rdpt_glyphs.data_source.data['position'] = originalRDPTData['position'].slice(0, newNumber)
                        rdpt_glyphs.data_source.data['score'] = originalRDPTData['score'].slice(0, newNumber)
                        rdpt_glyphs.data_source.data['strand'] = originalRDPTData['strand'].slice(0, newNumber)
                        rdpt_glyphs.data_source.change.emit()
                        """
        rdpt_javascript = CustomJS(
            args=dict(
                rdpt_glyphs=terminator_glyphs,
                rdpt_number=rdpt_number,
                source=ritterminator_dict,
            ),
            code=rdpt_javascript,
        )
        curdoc().js_on_event(DocumentReady, rdpt_javascript)
        fig.js_on_event("reset", rdpt_javascript)
        wigits.append((rdpt_number, rdpt_javascript))

        # Add table
        name = "Rho-Dependent Terminators"
        rdpt_table = generate_bokeh_table(rdpt, name)
        tables[name] = rdpt_table

    if ridpt:
        ridpt = sorted(ridpt, key=lambda x: x.conf, reverse=True)
        terminator_dict = {"x": [], "y": [], "position": [], "score": [], "strand": []}
        for terminator in ridpt:
            t_shape = (
                (0, 0),
                (20, 0),
                (20, 350),
                (60, 350),
                (60, 450),
                (-40, 450),
                (-40, 350),
                (0, 350),
                (0, 0),
            )
            scaled_t = [
                (
                    x[0] * len(cryptresult.sequence) / 10000 * GLYPH_SCALING,
                    x[1] * GLYPH_SCALING,
                )
                for x in t_shape
            ]
            t_shape = scaled_t
            scaled_position = terminator.start
            if terminator.strand == "-":
                t_shape = [(-x[0], -x[1]) for x in t_shape]
            xs = [scaled_position + x[0] for x in t_shape]
            if terminator.strand == "+":
                ys = [x[1] + forward_height for x in t_shape]
            else:
                ys = [x[1] + reverse_height for x in t_shape]
            terminator_dict["x"].append(xs)
            terminator_dict["y"].append(ys)
            terminator_dict["position"].append(f"{terminator.start}-{terminator.end}")
            terminator_dict["score"].append(terminator.conf)
            terminator_dict["strand"].append(terminator.strand)

        terminator_glyphs = fig.patches(
            "x",
            "y",
            color="blue",
            source=terminator_dict,
            alpha=0.5,
            line_color="black",
            line_width=1,
            y_range_name="y_range2",
        )
        terminator_glyphs_hover = HoverTool(
            renderers=[terminator_glyphs],
            tooltips=[
                ("Type", "Intrinsic Terminator"),
                ("Position", "@position"),
                ("Strand", "@strand"),
                ("Score", "@score"),
            ],
            visible=False,
        )
        fig.add_tools(terminator_glyphs_hover)

        # Add a widget to change the number of terminators shown
        ridpt_number = TextInput(
            title="Number of Intrinsic Terminators", value=str(number_of_glyphs)
        )
        ridpt_javascript = """
                        var originalRIDPTData = structuredClone(source)
                        var newNumber = Number(ridpt_number.value)
                        ridpt_glyphs.data_source.data['x'] = originalRIDPTData['x'].slice(0, newNumber)
                        ridpt_glyphs.data_source.data['y'] = originalRIDPTData['y'].slice(0, newNumber)
                        ridpt_glyphs.data_source.data['position'] = originalRIDPTData['position'].slice(0, newNumber)
                        ridpt_glyphs.data_source.data['score'] = originalRIDPTData['score'].slice(0, newNumber)
                        ridpt_glyphs.data_source.data['strand'] = originalRIDPTData['strand'].slice(0, newNumber)
                        ridpt_glyphs.data_source.change.emit()
                        """
        ridpt_javascript = CustomJS(
            args=dict(
                ridpt_glyphs=terminator_glyphs,
                ridpt_number=ridpt_number,
                source=terminator_dict,
            ),
            code=ridpt_javascript,
        )
        curdoc().js_on_event(DocumentReady, ridpt_javascript)
        fig.js_on_event("reset", ridpt_javascript)
        wigits.append((ridpt_number, ridpt_javascript))

        # Add table
        name = "Intrinsic Terminators"
        ridpt_table = generate_bokeh_table(ridpt, name)
        tables[name] = ridpt_table

    # Add the expressed CDSs
    expressed_CDSs = cryptresult.translation_sites
    color_bar_figure = None
    if expressed_CDSs:
        # Sort the expressed_CDS by the difference between the start and end
        def sort_algorythm(x):
            candidate_1 = x.end - x.start
            candidate_2 = x.start - x.end
            return max(candidate_1, candidate_2)

        expressed_CDSs = sorted(expressed_CDSs, key=sort_algorythm, reverse=True)

        if view_format == "mirrored":
            # Get the forward strand CDSs
            expressed_CDSs_fwd = [x for x in expressed_CDSs if x.strand == "+"]
            boxes_fwd, highest_tir_fwd = plot_boxes(expressed_CDSs_fwd)
            highest_y_pos = max(
                [(x[0] + x[1]) for x in zip(boxes_fwd["y"], boxes_fwd["h"])]
            )

            # Get the reverse strand CDSs
            expressed_CDSs_rev = [x for x in expressed_CDSs if x.strand == "-"]
            boxes_rev, highest_tir_rev = plot_boxes(expressed_CDSs_rev)
            highest_y_neg = max(
                [(x[0] + x[1]) for x in zip(boxes_rev["y"], boxes_rev["h"])]
            )
        elif view_format == "stacked":
            expressed_CDSs_fwd = [x for x in expressed_CDSs]
            boxes_fwd, highest_tir_fwd = plot_boxes(expressed_CDSs_fwd)
            highest_y_pos = max(
                [(x[0] + x[1]) for x in zip(boxes_fwd["y"], boxes_fwd["h"])]
            )
            boxes_rev, highest_tir_rev = (
                {
                    "x": [],
                    "y": [],
                    "w": [],
                    "h": [],
                    "rbs_strength": [],
                    "burden": [],
                    "position": [],
                    "strand": [],
                },
                0,
            )
            highest_y_neg = 0
            expressed_CDSs_rev = []
        else:
            raise ValueError(
                f'view_format must be "mirrored" or "stacked". {view_format} is not a valid option.'
            )

        highest_tir = max(highest_tir_fwd, highest_tir_rev)

        # Calculate the space needed for the annotations
        total_range = highest_y_pos + highest_y_neg
        room_on_graph_for_annotations = (total_range * (ANNOTATION_SPACE / 100)) / (
            1 - (ANNOTATION_SPACE / 100)
        )
        total_range = total_range + room_on_graph_for_annotations

        change_each_axis_by = 0

        # Flip the negative y values and make space for annotations
        for i in range(len(boxes_fwd["y"])):
            boxes_fwd["y"][i] = boxes_fwd["y"][i]
            boxes_fwd["y"][i] += change_each_axis_by
        for i in range(len(boxes_rev["y"])):
            boxes_rev["y"][i] = -boxes_rev["y"][i]
            boxes_rev["y"][i] -= change_each_axis_by

        # Plot the forward and the reverse (if they exist)
        cmap = linear_cmap(
            field_name="burden", palette=viridis(256), low=0, high=highest_tir
        )
        if forward_exists:
            source = ColumnDataSource(boxes_fwd)
            rectangles_fwd = fig.rect(
                x="x",
                y="y",
                width="w",
                height="h",
                source=source,
                color=cmap,
                line_color="black",
                line_width=1,
                y_range_name="y_range3",
                x_range_name="x_range3",
            )
            rectangles_fwd_hover = HoverTool(
                renderers=[rectangles_fwd],
                tooltips=[
                    ("Type", "Translation"),
                    ("Position", "@position"),
                    ("Strand", "@strand"),
                    ("RBS Strength", "@rbs_strength"),
                    ("Burden", "@burden"),
                ],
                visible=False,
            )
            fig.add_tools(rectangles_fwd_hover)
        if reverse_exists:
            source = ColumnDataSource(boxes_rev)
            rectangles_rev = fig.rect(
                x="x",
                y="y",
                width="w",
                height="h",
                source=source,
                color=cmap,
                line_color="black",
                line_width=1,
                y_range_name="y_range4",
                x_range_name="x_range3",
            )
            rectangles_rev_hover = HoverTool(
                renderers=[rectangles_rev],
                tooltips=[
                    ("Type", "Translation"),
                    ("Position", "@position"),
                    ("Strand", "@strand"),
                    ("RBS Strength", "@rbs_strength"),
                    ("Burden", "@burden"),
                ],
                visible=False,
            )
            fig.add_tools(rectangles_rev_hover)

        boxes = {
            key: boxes_fwd.get(key, 0) + boxes_rev.get(key, 0)
            for key in set(boxes_fwd) | set(boxes_rev)
        }

        source = ColumnDataSource(boxes)

        if forward_exists:
            color_bar = rectangles_fwd.construct_color_bar(
                padding=0, ticker=fig.xaxis.ticker, formatter=fig.xaxis.formatter
            )
        else:
            color_bar = rectangles_fwd.construct_color_bar(
                padding=0, ticker=fig.xaxis.ticker, formatter=fig.xaxis.formatter
            )

        color_bar_figure = figure(
            title="Burden",
            title_location="right",
            width=100,
            height=750,
            y_range=(0, highest_tir),
            toolbar_location=None,
            min_border=0,
            outline_line_color=None,
        )
        color_bar_figure.add_layout(color_bar, "right")
        color_bar_figure.title.align = "center"
        color_bar_figure.title.text_font_size = FONTSIZE
        color_bar_figure.title.text_font = FONT
        color_bar_figure.title.text_font_style = "bold"
        color_bar.major_label_text_font = FONT
        color_bar.major_label_text_font_size = FONTSIZE
        silence(warning=MISSING_RENDERERS)

        max_y_pos = TextInput(
            title="Max Y (Top track)", value=str(math.ceil(highest_y_pos))
        )
        max_y_neg = TextInput(
            title="Max Y (Bottom track)", value=str(math.ceil(highest_y_neg))
        )

        # fix the Y axis tickers
        ticker_locations = [n for n in range(0, 1000 * 1000, tick_frequency)]
        shownaxis.ticker = ticker_locations + [
            y * -1 - room_on_graph_for_annotations for y in ticker_locations
        ]
        ticker_labels = {n: str(n) for n in ticker_locations}
        ticker_labels = {
            **ticker_labels,
            **{
                y * -1 - room_on_graph_for_annotations: str(y) for y in ticker_locations
            },
        }
        # Round the labels to the nearest int
        fig.yaxis.major_label_overrides = {
            k: str(abs(int(float(v)))) for k, v in ticker_labels.items()
        }

        max_y_js = CustomJS(
            args=dict(
                plotRange1=fig.y_range,
                plotRangeAnnotations=fig.extra_y_ranges["y_range2"],
                plotRangeTop=fig.extra_y_ranges["y_range3"],
                plotRangeBot=fig.extra_y_ranges["y_range4"],
                max_y_top=max_y_pos,
                may_y_bot=max_y_neg,
                annotation_space=ANNOTATION_SPACE,
                annotation_depth=annotation_depth,
                figticker=shownaxis,
                tickFrequency=tick_frequency,
            ),
            code="""
    // Perform the math necessary to keep the annotation section consistant
    var newMaxYPos = Number(max_y_top.value)
    var newMaxYNeg = Number(may_y_bot.value)
    var total_range = newMaxYPos + newMaxYNeg
    var AnnotationSpace = Number(annotation_space)
    var AnnotationSpace = (total_range * (AnnotationSpace / 100))  / (1-(AnnotationSpace / 100))
    var total_space = newMaxYPos + newMaxYNeg + AnnotationSpace

    plotRangeTop.start = newMaxYPos-total_space
    plotRangeTop.end = newMaxYPos

    plotRangeBot.start = -newMaxYNeg
    plotRangeBot.end = total_space-newMaxYNeg

    var annotation_percent = AnnotationSpace / total_space
    var annotation_depth = Number(annotation_depth)
    var top_percent = newMaxYPos / total_space
    var bot_percent = newMaxYNeg / total_space

    var ptsPerPercent = -annotation_depth / annotation_percent
    var toppts = ptsPerPercent * top_percent
    var botpts = ptsPerPercent * bot_percent

    plotRangeAnnotations.start = -botpts+annotation_depth
    plotRangeAnnotations.end = toppts

    console.log(figticker)

    // Update the position of the reverse strand tickers
    let tickerLocations = Array.from({length: 1000000 / tickFrequency}, (_, n) => n * tickFrequency);
    const yaxisTicker = Array(...tickerLocations, ...tickerLocations.map(y => -y - AnnotationSpace));
    figticker.ticker.ticks = yaxisTicker

    // Update the labels of the reverse strand tickers
    const ticker_labels = new Map();
    for (let n of tickerLocations) {
        ticker_labels.set(n, String(Math.round(n)));
        ticker_labels.set(-n - AnnotationSpace, String(Math.round(n)));
    }
    figticker.major_label_overrides = ticker_labels

""",
        )
        curdoc().js_on_event(DocumentReady, max_y_js)
        fig.js_on_event("reset", max_y_js)

        wigits.append((max_y_pos, max_y_js))
        wigits.append((max_y_neg, max_y_js))

        max_burden = TextInput(title="Max burden", value=str(math.ceil(highest_tir)))
        max_burden_js = CustomJS(
            args=dict(color_bar=color_bar),
            code="""
    var newMax = Number(cb_obj.value)
    color_bar.color_mapper.high = newMax
""",
        )
        wigits.append((max_burden, max_burden_js))

        # Add table
        expressed_CDSs = sorted(expressed_CDSs, key=lambda x: x.burden, reverse=True)
        name = "Expressed CDSs"
        expressed_CDSs_table = generate_bokeh_table(expressed_CDSs, name)
        tables[name] = expressed_CDSs_table

        if not show_small:
            # Javascript to hide CDSs less than 45 bases
            js_string = """
        var newMin = 45
        rectangles.visible = true
        for (var i = 0; i < rectangles.data_source.data['w'].length; i++) {
            if (rectangles.data_source.data['w'][i] < newMin) {
                rectangles.data_source.data['w'][i] = 0
                rectangles.data_source.data['h'][i] = 0
            }
        }
        rectangles.data_source.change.emit()
    """
            if forward_exists:
                hide_short_CDSs_js_fwd = CustomJS(
                    args=dict(rectangles=rectangles_fwd), code=js_string
                )
                curdoc().on_event(DocumentReady, hide_short_CDSs_js_fwd)
                fig.js_on_event("reset", hide_short_CDSs_js_fwd)
            if reverse_exists:
                hide_short_CDSs_js_rev = CustomJS(
                    args=dict(rectangles=rectangles_rev), code=js_string
                )
                curdoc().on_event(DocumentReady, hide_short_CDSs_js_rev)
                fig.js_on_event("reset", hide_short_CDSs_js_rev)

    # Build a DIV above the plot that contains the name of the plot and the total burden
    if cryptresult.name:
        name_div = Div(
            text=f'<h1 style="font-family: {FONT}; font-size: {FONTSIZE}">{cryptresult.name}</h1>'
        )
    else:
        name_div = Div()
    burden_div = Div(
        text=f'<h2 style="font-family: {FONT}; font-size: {FONTSIZE}">Total Burden: {round(cryptresult.burden, 2)}</h2>'
    )

    widgets_to_add = []
    for wigit, js in wigits:
        wigit.js_on_change("value", js)
        widgets_to_add.append(wigit)

    # Specify we want the svg backend
    fig.output_backend = "svg"

    # Add the plot to the document and add a title
    widgets = row(*widgets_to_add, styles={"margin": "0 auto", "align-items": "center"})
    if color_bar_figure:
        layout = row(fig, color_bar_figure)
    else:
        layout = fig

    # layout = column(layout, widgets, styles={'margin': '0 auto', 'align-items': 'center'})
    layout = column(
        name_div,
        burden_div,
        layout,
        widgets,
        styles={"margin": "0 auto", "align-items": "center"},
    )
    layout.sizing_mode = "scale_width"

    # Add the tables
    if tables:
        if len(tables) > 1:
            # Make buttons that show the appropriat tables
            buttons = []
            button_label = Div(text="Tables:")
            for table_name, table_value in tables.items():
                button = bokeh.models.widgets.Button(label=table_name)
                button.js_on_click(
                    CustomJS(
                        args=dict(tables=tables, label=table_name),
                        code="""
                tables[label].visible = true;
                for (var key in tables) {
                    if (key != label){
                        tables[key].visible = false;
                    }
                }
                """,
                    )
                )
                buttons.append(button)
                table_value.visible = False
            buttons = row(
                button_label,
                *buttons,
                styles={"margin": "0 auto", "align-items": "center"},
            )
            layout = column(
                layout, buttons, styles={"margin": "0 auto", "align-items": "center"}
            )
            layout = column(
                layout,
                *tables.values(),
                styles={"margin": "0 auto", "align-items": "center"},
            )
        else:
            layout = column(
                layout, *tables, styles={"margin": "0 auto", "align-items": "center"}
            )

    if filename:
        script, fig = components(layout)
        export_html(script, fig, filename + "_graph.html")
        # copy the assets folder to the output folder if it doesn't exist
        if os.path.exists(os.path.join(os.path.dirname(filename), "assets")):
            shutil.rmtree(os.path.join(os.path.dirname(filename), "assets"))
        shutil.copytree(
            os.path.join(os.path.dirname(__file__), "assets"),
            os.path.join(os.path.dirname(filename), "assets"),
        )
        for widget, js in wigits:
            widget.visible = False
        return filename

    return fig


def generate_bokeh_table(datalist, name) -> DataTable:
    # Generate a bokeh table from a list of named tuples
    column_names = datalist[0]._fields
    table_name = datalist[0].__class__.__name__
    data = {column_name: [] for column_name in column_names}
    for datarow in datalist:
        for column_name in column_names:
            data[column_name].append(getattr(datarow, column_name))
    source = ColumnDataSource(data)

    columns = []
    for column_name in column_names:
        if column_name in ["score", "c_over_g", "tail_score", "burden", "dG"]:
            formatted_column = TableColumn(
                field=column_name,
                title=column_name,
                formatter=NumberFormatter(format="‘0,0.0’", rounding="round"),
            )
        else:
            formatted_column = TableColumn(field=column_name, title=column_name)
        columns.append(formatted_column)

    name_div = Div(text=f"<h1>{name}</h1>")
    table = DataTable(
        source=source, columns=columns, name=table_name, width=1500, editable=True
    )
    table = column(name_div, table)
    return table


def export_html(script, div, filename):
    bokeh_version = bokeh.__version__
    template_path = os.path.join(
        os.path.dirname(__file__), "assets", "html_template.html"
    )
    asset_location = os.path.join(os.path.dirname(__file__), "assets")
    with open(template_path, "r", encoding="utf-8") as f:
        template = f.read()
    template = template.replace("{{bokehscript}}", script)
    template = template.replace("{{bokehdiv}}", div)
    template = template.replace("{{bokehversion}}", bokeh_version)
    template = template.replace("{{assetlocation}}", asset_location)
    with open(filename, "w", encoding="utf-8") as f:
        f.write(template)
    return filename


# EOF
