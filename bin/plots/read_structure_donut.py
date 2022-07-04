import imp
from subprocess import check_output
from collections import Counter
import os
from math import pi

import pandas as pd
import pysam

# from plotly.subplots import make_subplots
# import plotly.graph_objects as go
# import plotly.express as px

from bokeh.io import show, save, output_file
from bokeh.plotting import figure
from bokeh.transform import cumsum
from bokeh.layouts import row, column
from bokeh.models import Title, Div, ColumnDataSource


def _bam_to_df(bam, chr=None, start=None, stop=None):
    """
    Convert a bam into a pandas dataframe for the columns we need.
    """
    qname = []
    flag = []
    rname = []
    pos = []
    mapq = []
    cigar = []
    original_length = []

    for read in bam.fetch(chr, start, stop):
        qname.append(read.qname)
        flag.append(read.flag)
        rname.append(read.reference_name)
        pos.append(read.pos)
        mapq.append(read.mapq)
        cigar.append(read.cigarstring)
        original_length.append(read.infer_read_length())

    return pd.DataFrame(
        {
            "QNAME": qname,
            "FLAG": flag,
            "RNAME": rname,
            "POS": pos,
            "MAPQ": mapq,
            "CIGAR": cigar,
            "original_length": original_length,
        }
    )



def create_split_bam_table(bam_path):
    """
    Create a dataframe and set the properties correctly.
    """
    pysam.index(bam_path)
    aln_file = pysam.AlignmentFile(bam_path, "rb")
    split_table = _bam_to_df(aln_file)
    split_table["FLAG"] = split_table.FLAG.astype(int)
    split_table["POS"] = split_table.POS.astype(int)
    split_table["MAPQ"] = split_table.MAPQ.astype(int)
    split_table["original_length"] = split_table.original_length.astype(int)
    return split_table


def count_reads_and_bases(split_table):
    """
    Determine the type of read based on the present alignments.
    """

    def initialize_counters():
        concat_types = {
            "BB-I": 0,
            "BB-only": 0,
            "I-only": 0,
            "Unmapped": 0,
            "mBB-only": 0,  # multiple BB, no I
            "mI-only": 0,  # multiple I, no BB
            "BB-mI": 0,  # single BB, multiple I
            "mBB-I": 0,  # multiple BB, single I
            "mBB-mI": 0,  # multiple BB, multiple I
            "Unknown": 0,  # anything else
        }
        concat_type_stats = Counter(concat_types)
        concat_type_stats_by_bases = Counter(concat_types)
        return concat_type_stats, concat_type_stats_by_bases

    concat_type_stats, concat_type_stats_by_bases = initialize_counters()

    for qname, g in split_table.groupby("QNAME"):
        mapped = list(g.RNAME.unique())
        mapped_bb = [i for i in mapped if i.startswith("BB")]
        mapped_ins = [i for i in mapped if i not in mapped_bb]
        raw_read_length = g.original_length.iloc[0]

        # single element mapped
        if len(mapped) == 1:
            if mapped == ["*"]:
                concat_type_stats.update({"Unmapped"})
                concat_type_stats_by_bases.update({"Unmapped": raw_read_length})
            elif mapped[0] in mapped_bb:
                concat_type_stats.update({"BB-only"})
                concat_type_stats_by_bases.update({"BB-only": raw_read_length})
            else:
                concat_type_stats.update({"I-only"})
                concat_type_stats_by_bases.update({"I-only": raw_read_length})

        # two elements mapped
        elif len(mapped) == 2:
            if mapped_bb and mapped_ins:
                concat_type_stats.update({"BB-I"})
                concat_type_stats_by_bases.update({"BB-I": raw_read_length})
            elif mapped_bb and not mapped_ins:
                concat_type_stats.update({"BB-only"})
                concat_type_stats_by_bases.update({"BB-only": raw_read_length})
            else:
                concat_type_stats.update({"I-only"})
                concat_type_stats_by_bases.update({"I-only": raw_read_length})

        # complex mapping
        else:
            if not mapped_ins and len(mapped_bb) > 1:
                concat_type_stats.update({"mBB-only"})
                concat_type_stats_by_bases.update({"mBB-only": raw_read_length})

            elif not mapped_bb and len(mapped_ins) > 1:
                concat_type_stats.update({"mI-only"})
                concat_type_stats_by_bases.update({"mI-only": raw_read_length})

            elif len(mapped_ins) > 1 and len(mapped_bb) == 1:
                concat_type_stats.update({"BB-mI"})
                concat_type_stats_by_bases.update({"BB-mI": raw_read_length})

            elif len(mapped_bb) > 1 and len(mapped_ins) == 1:
                concat_type_stats.update({"mBB-I"})
                concat_type_stats_by_bases.update({"mBB-I": raw_read_length})

            elif len(mapped_bb) > 1 and len(mapped_ins) > 1:
                concat_type_stats.update({"mBB-mI"})
                concat_type_stats_by_bases.update({"mBB-mI": raw_read_length})
            else:
                concat_type_stats.update({"Unknown"})
                concat_type_stats_by_bases.update({"Unknown": raw_read_length})

    print(concat_type_stats)
    print(concat_type_stats_by_bases)
    return concat_type_stats, concat_type_stats_by_bases


def create_assembly_count_plot(split_table, output_file_name, my_title):
    """
    Create a barplot for the assembly statistics.
    """
    ##Segments count
    _d = Counter(split_table.RNAME.values)
    segments = (
        pd.DataFrame.from_dict(_d, orient="index")
        .reset_index()
        .rename(columns={"index": "CHROM", 0: "count"})
    )
    print(segments)
    source = ColumnDataSource(data=segments)

    output_file(filename=output_file_name, title=my_title)
    p = figure(plot_height=500, plot_width=1000, x_range=segments.CHROM, title=my_title)
    p.vbar(x="CHROM", top="count", width=0.8, source=source)
    p.xaxis.major_label_orientation = "vertical"
    save(p)


def create_readtype_donuts(concat_type_stats, concat_type_stats_by_bases, plot_file, plot_title):
    """
    Create a donut plot for the read structures.
    """

    concat_type_colors = {
        "BB-I": "DodgerBlue",  # perfect
        "BB-only": "Crimson",  # waste
        "I-only": "MediumSlateBlue",  # informing but missing barcode
        "Unmapped": "Gray",  # waste
        "mBB-only": "Tomato",  # waste
        "mI-only": "DarkOrchid",  # informing but missing barcode
        "BB-mI": "RoyalBlue",  # OK
        "mBB-I": "Navi",  # OK
        "mBB-mI": "SkyBlue",  # OK
        "Unknown": "Gold",  # Waste
    }

    def _calculate_angle_and_color(stats, concat_type_colors):
        _df1 = pd.DataFrame.from_dict(stats, orient="index").reset_index()
        _df1.columns = ["type", "count"]
        _df1["angle"] = _df1["count"] / _df1["count"].sum() * 2 * pi
        _df1["color"] = _df1.type.map(lambda x: concat_type_colors[x])

        return _df1


    def _plot_donut(
        data,
        donut_plot_height,
        donut_plot_width,
        donut_plot_x_range,
        donut_plot_y_range,
        subtitle
    ):
        # data = ColumnDataSource(data)

        # create empty plot object
        donut_plot = figure(
            plot_height=donut_plot_height,
            plot_width=donut_plot_width,
            title=subtitle,
            tools="hover",
            tooltips="@type: @count",
            x_range=donut_plot_x_range,
            y_range=donut_plot_y_range,
            toolbar_location=None,
        )
        # add wedges per type
        donut_plot.annular_wedge(
            x=0,
            y=1,
            inner_radius=0.2,
            outer_radius=0.45,
            direction="anticlock",
            start_angle=cumsum("angle", include_zero=True),
            end_angle=cumsum("angle"),
            line_color="white",
            fill_color="color",
            legend_group="type",
            source=data,
        )

        # remove chart elements
        donut_plot.axis.axis_label = None
        donut_plot.axis.visible = False
        donut_plot.grid.grid_line_color = None
        return donut_plot


    output_file(filename=plot_file, title=plot_title)
    _df1 = _calculate_angle_and_color(concat_type_stats, concat_type_colors)
    _df2 = _calculate_angle_and_color(concat_type_stats_by_bases, concat_type_colors)

    donut_plot_height = 400
    donut_plot_width = 600
    donut_plot_x_range = (-0.6, 1.4)
    donut_plot_y_range = (0, 2)
    p1 = _plot_donut(
        _df1,
        donut_plot_height,
        donut_plot_width,
        donut_plot_x_range,
        donut_plot_y_range,
        "per read"
    )
    p2 = _plot_donut(
        _df2,
        donut_plot_height,
        donut_plot_width,
        donut_plot_x_range,
        donut_plot_y_range,
        "per base"
    )

    donut_row = row(p1, p2)
    donut_plot = column(Div(text=f"<h1>{plot_title}</h1>"), donut_row)
    save(donut_plot)



def main(bam, assembly_plot_file, assembly_plot_title, donut_plot_file, donut_plot_title):
    """
    Extract data from bam and plot it using two functions that further select data
    """
    df = create_split_bam_table(bam)
    read_counts, base_counts = count_reads_and_bases(df)
    create_assembly_count_plot(df, assembly_plot_file, assembly_plot_title)
    create_readtype_donuts(read_counts, base_counts, donut_plot_file, donut_plot_title)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Create plots for a split read bam."
    )

    parser.add_argument("bam_file")
    parser.add_argument("assembly_plot")
    parser.add_argument("donut_plot_readstructure")
    args = parser.parse_args()


    # split_bam = "/media/dami/a2bc89fb-be6b-4e23-912a-0c7137cd69ad/results/Cyclomics/000014/CycasConsensus/Minimap2AlignAdaptiveParameterized/FAT55621_pass_61af019d_23_filtered.bam"
    # output_file_name = "tmp.html"
    # title = "my title"
    # assert os.path.exists(split_bam)
    # main(split_bam, "assem.html", 'assem count', "Donuts.html", 'Donut plot 123')

    main(args.bam_file, args.assembly_plot, "Segment count per assembly in reference", args.donut_plot_readstructure, "Read structure by both read and base count")
