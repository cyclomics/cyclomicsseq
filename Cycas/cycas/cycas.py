import sys

import click
from loguru import logger
from src import app

logger.remove()
logger.add(sys.stderr, level="INFO")

old = False

# Command Group, all others could initiate from this
@click.group(name="tools")
def cli():
    """Tool related commands"""


@cli.command(
    name="consensus", help="full consensus generation, creates consensus Fastq's."
)
@click.option(
    "--bam-file",
    required=True,
    help="path to the bamfile containing the secondary alignment files",
    type=click.Path(exists=True),
)
@click.option(
    "--output",
    required=True,
    help="fastq file to write the output to",
    type=click.Path(exists=False),
)
@click.option(
    "--metadata-json",
    help="Path to write metadata in json format to.",
    required=False,
    default=None,
    type=click.Path(exists=False),
)
@click.option("--plot-readtypes/--no-plot-readtypes", default=False)
@click.option("--plot-piechart/--no-plot-piechart", default=False)
@click.option("--show-class-counts/--no-show-class-counts", default=False)
@click.option(
    "--limit-calls",
    help="Limit the amount of reads to classify, usefull for testing. Set to 0 to disable.",
    default=0,
)
def full_consensus(
    bam_file,
    output,
    metadata_json,
    plot_readtypes,
    plot_piechart,
    show_class_counts,
    limit_calls,
):
    """
    Aim:
    Given a Samfile, classify and create consensus for each read and if possible add the barcode.
    Using:
        1. Directional chromosomes
        2. starting site within chromosome
        3. continuity of the alignment

    Additional:
        1. Output fastq of partial alignments
        2. Create spagettograms
            - per type
            - specific read
        3. piechart of total classes
        4. use set of filters
        5. option to change the consensus caller based on class

    Output:
        1. consensus Fastq with barcode in name
        2. metadata json
        3. Optional: fragmented alignments
    """
    click.echo("Cycas: Full consensus generation")

    app.create_classifications_consensus(
        bam_file=bam_file,
        output=output,
        plot_readtypes=plot_readtypes,
        plot_piechart=plot_piechart,
        metadata_json=metadata_json,
        limit_calls=limit_calls,
    )


@cli.command(
    name="single-read", help="classify and generate consensus for a single read"
)
@click.option(
    "--bam-file",
    required=True,
    help="path to the bamfile containing the secondary alignment files.",
    type=click.Path(exists=True),
)
@click.option(
    "--output",
    required=True,
    help="fastq file to write the output to",
    type=click.Path(exists=False),
)
@click.option("--read-name", help="Name of the read in the bam file.", required=True)
@click.option(
    "--metadata-json",
    help="Path to write metadata in json format to.",
    required=False,
    default=None,
    type=click.Path(exists=False),
)
def single_read(bam_file, output, read_name, metadata_json):
    click.echo("Cycas: Classification and splitting a single read.")
    app.classify_single_read(bam_file, output, read_name, metadata_json)


@cli.command(
    name="split",
    help="classify and separate the reads, keeps them aligned no consensus generation.",
)
@click.option(
    "--bam-file",
    required=True,
    help="path to the bamfile containing the secondary alignment files",
    type=click.Path(exists=True),
)
@click.option(
    "--output",
    required=True,
    help="fastq file to write the output to",
    type=click.Path(exists=False),
)
@click.option(
    "--metadata-json",
    help="Path to write metadata in json format to.",
    required=False,
    default=None,
    type=click.Path(exists=False),
)
@click.option("--plot-readtypes/--no-plot-readtypes", default=False)
@click.option("--plot-piechart/--no-plot-piechart", default=False)
@click.option(
    "--limit-calls",
    help="Limit the amount of reads to classify, usefull for testing. Set to 0 to disable.",
    default=0,
)
def split_reads(
    bam_file,
    output,
    plot_readtypes,
    plot_piechart,
    show_class_counts,
    create_classification_detail_json,
    create_metadata_json,
    limit_calls,
):
    click.echo("Cycas: Classification and splitting.")
    raise NotImplementedError(
        "splitting is currently not implemented in this version of Cycas."
    )
    app.separate_reads(
        bam_file=bam_file,
        output=output,
        plot_readtypes=plot_readtypes,
        plot_piechart=plot_piechart,
        create_classification_detail_json=create_classification_detail_json,
        create_metadata_json=create_metadata_json,
        limit_calls=limit_calls,
    )


if __name__ == "__main__":
    cli()
