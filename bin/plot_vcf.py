import io
import pandas as pd


from bokeh.io import save, output_file
from bokeh.plotting import figure, show
from bokeh.layouts import row, column
from bokeh.models import HoverTool


def read_vcf(path):
    with open(path, "r") as f:
        lines = [l for l in f if not l.startswith("##")]
    df = pd.read_csv(
        io.StringIO("".join(lines)),
        dtype={
            "#CHROM": str,
            "POS": int,
            "ID": str,
            "REF": str,
            "ALT": str,
            "QUAL": str,
            "FILTER": str,
            "INFO": str,
        },
        sep="\t",
    ).rename(columns={"#CHROM": "CHROM"})

    formats = df.FORMAT[0].split(":")
    for i, fmt in enumerate(formats):
        df[fmt] = df.Sample1.apply(
            lambda x: float((x.split(":")[i] if (x.split(":")[i]) else 0))
        )
    del df["FORMAT"]
    del df["Sample1"]
    return df


def update_vcf_with_truth(true_vcf, df):
    true_df = read_vcf(true_vcf)
    true_df = true_df[["CHROM", "POS"]]
    true_df["Variant"] = 1
    df = pd.merge(
        df, true_df, "left", left_on=["CHROM", "POS"], right_on=["CHROM", "POS"]
    )
    df["Variant"].fillna(0, inplace=True)

    return df


def make_scatter_plots(data):
    hover = HoverTool(
        tooltips=[
            ("position", "@CHROM : @POS"),
            ("VAF", "@VAF"),
            ("Frequency", "@FREQ"),
            ("forward ratio", "@FWDR"),
            ("reverse ratio", "@REVR"),
            ("Same base", "@SAME"),
            ("Depth", "@DP"),
            ("Depth after Qfilter", "@DPQ"),
        ]
    )

    # Pos scatter vaf
    p_vaf = figure(title="Positional VAF")
    p_vaf.xaxis.axis_label = "Position"
    p_vaf.yaxis.axis_label = "VAF ratio"
    p_vaf.scatter("POS", "VAF", source=data)
    p_vaf.add_tools(hover)
    # show(p_vaf)

    # fwd ratio vs reverse ratio
    p_ratio = figure(title="Support ratio's per direction")
    p_ratio.xaxis.axis_label = "forward  ratio"
    p_ratio.yaxis.axis_label = "reverse ratio"
    p_ratio.scatter("FWDR", "REVR", source=data)
    p_ratio.add_tools(hover)
    # show(p)

    # Pos scatter FREQ
    p_freq = figure(title="Positional FREQ")
    p_freq.xaxis.axis_label = "Position"
    p_freq.yaxis.axis_label = "FREQ"
    p_freq.scatter("POS", "FREQ", source=data)
    p_freq.add_tools(hover)
    # show(p)

    # Pos scatter FREQ
    p_depth = figure(title="Depth pre and post base-quality filtering")
    p_depth.xaxis.axis_label = "Position"
    p_depth.yaxis.axis_label = "Depth"
    p_depth.scatter("POS", "DP", source=data, color="orange", legend_label="depth")
    p_depth.scatter("POS", "DPQ", source=data, legend_label="depth after Q filtering")
    p_depth.add_tools(hover)
    p_depth.legend.location = "bottom_center"

    return column([row([p_vaf, p_ratio]), row([p_depth, p_freq])])


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Create hist plot from a regex for fastq and fastq.gz files."
    )

    parser.add_argument("vcf_file")
    parser.add_argument("plot_file")
    args = parser.parse_args()
    # vcf = '/home/dami/projects/variantcalling/depth/datasets/000010_v4_bb41/2000/potential_vars_1000.vcf'
    # true_vcf = '/home/dami/projects/variantcalling/depth/datasets_cyclomics/000010_v4_bb41/Native/real_variants.vcf'

    data = read_vcf(args.vcf_file)
    # print(data)
    plot = make_scatter_plots(data)
    output_file(args.plot_file, title="variant plots")
    save(plot)
