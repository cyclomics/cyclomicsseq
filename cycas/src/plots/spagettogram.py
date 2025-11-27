import logging
import re
import time
from collections import OrderedDict, defaultdict
from pathlib import Path

from bokeh.io import output_file, save
from bokeh.layouts import column
from bokeh.palettes import Category10
from bokeh.plotting import figure, show
from loguru import logger


def process_forward_cigar(fig, cigars, pos_y):
    pos_x = 0
    aln_start = 0
    found_start = False
    y_vals = []

    for aln_type, amount in cigars:
        amount = int(amount)
        y_line = (pos_y, pos_y + amount)
        x_line = (pos_x, pos_x + amount)

        if aln_type == "M":
            if not found_start:
                found_start = True
                aln_start = pos_x

            fig.line(x_line, y_line, line_width=3, color="orange")
            pos_y += amount
            pos_x += amount

        elif aln_type == "I":
            # extra read nucleotides, so no y progress
            y_line = (pos_y, pos_y)
            logging.debug(f"drawing line for {x_line} ,{y_line}")
            fig.line(x_line, y_line, line_width=20, color="green")
            pos_x += amount

        elif aln_type == "D":
            # missing read nucleotides, so no x progress
            x_line = (pos_x, pos_x)
            fig.line(x_line, y_line, line_width=20, color="red")
            pos_y += amount

        elif aln_type in ["S", "H"]:
            pos_x += amount
        y_vals.append(pos_y)

    return aln_start, min(y_vals), max(y_vals)


def process_reverse_cigar(fig, cigars, pos_y):
    # slope down, progress X
    pos_x = 0
    aln_start = 0
    found_start = False
    y_vals = [pos_y]

    # Minimap needs its cigar strings reversed for
    position_adjustment_y = [x for x in cigars if x[0] not in ["S", "H"]]
    position_adjustment_y = [x[1] for x in position_adjustment_y]
    position_adjustment_y = sum(position_adjustment_y)
    pos_y += position_adjustment_y

    for aln_type, amount in cigars:
        y_line = (pos_y - amount, pos_y)
        x_line = (pos_x + amount, pos_x)

        if aln_type in ["S", "H"]:
            pos_x += amount
            if not found_start:
                found_start = True
                aln_start = pos_x

        elif aln_type == "M":
            fig.line(x_line, y_line, line_width=3, color="skyblue")
            pos_y -= amount
            pos_x += amount

        elif aln_type == "I":
            # extra read nucleotides, so no y progress
            y_line = (pos_y, pos_y)
            fig.line(x_line, y_line, line_width=20, color="green")
            pos_x += amount

        elif aln_type == "D":
            # missing read nucleotides, so no x progress
            x_line = (pos_x, pos_x)
            fig.line(x_line, y_line, line_width=20, color="red")
            pos_y -= amount

        y_vals.append(pos_y)

    return aln_start, min(y_vals), max(y_vals)


def update_chrom_minmax(chrom, min_val, max_val):
    """
    update a dict that contains min and max values observed per chromosome
    """
    if len(chrom) > 0:
        new_min = min((chrom[0], min_val))
        new_max = max((chrom[1], max_val))
        return (new_min, new_max)
    else:
        return (min_val, max_val)


def create_new_figure(plots, chrom):
    """
    helper function for new figures with shared x axis
    """
    figure_options = dict(width=1000, height=350, title=chrom)

    if len(plots.keys()) > 0:
        for key, value in plots.items():
            leading_x_axis = value.x_range
            break

        fig = figure(**figure_options, x_range=leading_x_axis)
    else:
        fig = figure(**figure_options)

    fig.xaxis.axis_label = "Nucleotide position in raw read"
    fig.yaxis.axis_label = "Nucleotide position in reference"
    fig.yaxis.formatter.use_scientific = False

    return fig


def create_spagettogram(alignments, title=None):
    """
    plot group with read position on the x axis and the reference position on the y axis
    """
    logger.info("\tCreating spagettogram")
    plots = OrderedDict()
    chroms_plotted = []
    x_coords = []

    # params needed for start lines
    chrom_minmax = defaultdict(list)
    aln_starts = []

    if len(alignments) == 0:
        logger.warning(f"Empty alignment list given to spagettogram plotter")
        return

    for aln in alignments:
        # we use the chromosome all over the place
        chrom = aln.alignment_chromosome

        if chrom not in plots.keys():
            plots[chrom] = create_new_figure(plots, chrom)

        fig = plots[chrom]
        # zero based exclusive, so we add one
        pos_y = aln.alignment_chromosome_start + 1

        # decide if we plot forwards or backwards
        if aln.alignment_direction == "+":
            processor = process_forward_cigar
        elif aln.alignment_direction == "-":
            processor = process_reverse_cigar

        start, y_min, y_max = processor(fig, aln.cigars, pos_y)
        chrom_minmax[chrom] = update_chrom_minmax(chrom_minmax[chrom], y_min, y_max)
        aln_starts.append(start)

    # Add the horizontal start lines to each chromosome plot
    for chrom, plt in plots.items():
        for start in aln_starts:
            logger.debug(f"Draw startsite at {start}")
            plt.line((start, start), chrom_minmax[chrom], color="grey")
        # add start and end lines
        plt.line(
            (alignments[0].raw_read_length, alignments[0].raw_read_length),
            chrom_minmax[chrom],
            color="black",
        )
        plt.line((0, 0), chrom_minmax[chrom], color="black")

    timestamp = int(time.time())
    readname = alignments[0].readname

    defined_tile = title if not None else readname
    if title:
        spagettogram_filename = f"plots/{timestamp}_spagettogram_{title}.html"
    else:
        spagettogram_filename = f"plots/{timestamp}_spagettogram_{readname}.html"

    Path(spagettogram_filename).parent.mkdir(exist_ok=True)

    output_file(filename=spagettogram_filename, title=defined_tile)
    save(column(list(plots.values())))
    logger.info(f"Created plot at {spagettogram_filename}")
