import os
import shutil
import sys

import click
from loguru import logger
from src import app
from src.config import ReadClassificationConfig

logger.remove()
logger.add(sys.stderr, level="INFO")


# Command Group, all others could initiate from this
@click.group(name="tools")
def cli():
    """Tool related commands"""


@cli.command(
    name="quick-consensus", help="full consensus generation, creates consensus Fastq's."
)
@click.option(
    "--input-bam",
    required=True,
    help="path to the bamfile containing the secondary alignment files",
    type=click.Path(exists=True),
)
@click.option(
    "--output-fastq",
    required=True,
    help="fastq file to write the output to",
    type=click.Path(exists=False),
)
@click.option(
    "--output-json",
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
@click.option(
    "--hdf5-path",
    help="Create custom hdf5 file for these consensus results.",
    default=None,
    type=click.Path(exists=False),
)
def cli_quick_consensus(
    input_bam,
    output_fastq,
    output_json,
    plot_readtypes,
    plot_piechart,
    show_class_counts,
    limit_calls,
    hdf5_path,
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

    app.create_quick_consensus(
        bam_file=input_bam,
        output=output_fastq,
        plot_readtypes=plot_readtypes,
        plot_piechart=plot_piechart,
        metadata_json=output_json,
        limit_calls=limit_calls,
        hdf5_file=hdf5_path,
    )


@cli.command(
    name="quick-consensus-sr", help="classify and generate consensus for a single read"
)
@click.option(
    "--input-bam",
    required=True,
    help="path to the bamfile containing the secondary alignment files.",
    type=click.Path(exists=True),
)
@click.option(
    "--output-fastq",
    required=True,
    help="fastq file to write the output to",
    type=click.Path(exists=False),
)
@click.option("--read-name", help="Name of the read in the bam file.", required=True)
@click.option(
    "--output-json",
    help="Path to write metadata in json format to.",
    required=False,
    default=None,
    type=click.Path(exists=False),
)
def cli_quick_consensus_sr(input_bam, output_fastq, read_name, output_json):
    click.echo("Cycas: Classification and splitting a single read.")
    app.create_quick_consensus_sr(input_bam, output_fastq, read_name, output_json)


@cli.command(
    name="classify",
    help="Split and classify reads from a BAM file using structural rules.",
)
@click.option(
    "--input-bam",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Path to a BAM file containing reads and alignments.",
)
@click.option(
    "--output-bam",
    required=True,
    default=None,
    type=click.Path(writable=True, dir_okay=False),
    help="Path to the output BAM file after read splitting.",
)
@click.option(
    "--output-json",
    required=True,
    default=None,
    type=click.Path(writable=True, dir_okay=False),
    help="Optional path to write read metadata in JSON format.",
)
@click.option(
    "--output-h5",
    required=False,
    default=None,
    type=click.Path(writable=True, dir_okay=False),
    help="Path to write the read MSA in HDF5 format.",
)
@click.option(
    "--limit-calls",
    default=0,
    type=int,
    show_default=True,
    help="Limit the number of reads to classify (0 disables the limit).",
)
@click.option(
    "--max-orientation-error-1d",
    default=0.1,
    type=float,
    show_default=True,
    help="Maximum orientation error allowed for 1D classification.",
)
@click.option(
    "--max-orientation-error-2d",
    default=0.2,
    type=float,
    show_default=True,
    help="Maximum orientation error allowed for 2D classification.",
)
@click.option(
    "--min-mapped-inserts",
    default=3,
    type=int,
    show_default=True,
    help="Minimum number of mapped inserts required for classification.",
)
@click.option(
    "--backbone-prefix",
    default="BB",
    type=str,
    show_default=True,
    help="Prefix for backbone reference names to ignore during classification.",
)
@click.option(
    "--consecutiveness-method",
    default="new_occurrence",
    type=click.Choice(
        [
            "new_occurrence",
        ],
        case_sensitive=False,
    ),
    show_default=True,
    help="Method to calculate whether a segment appears consecutively in the original concatemer.",
)
def cli_classify(
    input_bam,
    output_bam,
    output_json,
    output_h5,
    limit_calls,
    max_orientation_error_1d,
    max_orientation_error_2d,
    min_mapped_inserts,
    backbone_prefix,
    consecutiveness_method,
):
    """
    CLI entry point for Cycas structural read classification.

    Processes a BAM file and applies structural classification rules to the reads.
    Outputs metadata optionally as JSON, and prints classification statistics.

    Args:
        bam_file: Path to a BAM file to classify.
        output_bam: Path to the output BAM file after read splitting.
        output_json: Optional path to write metadata.
        limit_calls: Limit number of reads processed (0 = no limit).
        max_orientation_error_1d: Max 1D orientation error allowed.
        max_orientation_error_2d: Max 2D orientation error allowed.
        min_mapped_inserts: Minimum inserts for classification.
        backbone_prefix: Prefix used to exclude backbone references.
        consecutiveness_method: Method to calculate whether a segment appears consecutively in the original concatemer.
    """
    click.echo("Cycas: 3D Classification and splitting.")

    config = ReadClassificationConfig(
        max_orientation_error_1D=max_orientation_error_1d,
        max_orientation_error_2D=max_orientation_error_2d,
        min_mapped_inserts=min_mapped_inserts,
        backbone_prefix=backbone_prefix,
        consecutiveness_method=consecutiveness_method,
    )

    app.create_classifications(
        bam_file=input_bam,
        output_bam=output_bam,
        metadata_json=output_json,
        hdf5_file=output_h5,
        limit_calls=limit_calls,
        config=config,
    )


@cli.command(
    name="consensus",
    help="Split, classify, call consensus and calibrate concatemeric reads.",
)
@click.option(
    "--input-bam",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Path to a BAM file containing reads and alignments.",
)
@click.option(
    "--output-fastq",
    required=True,
    default=None,
    type=click.Path(writable=True, dir_okay=False),
    help="Path to the output FASTQ file after read splitting.",
)
@click.option(
    "--output-json",
    required=True,
    default=None,
    type=click.Path(writable=True, dir_okay=False),
    help="Optional path to write read metadata in JSON format.",
)
@click.option(
    "--limit-calls",
    default=0,
    type=int,
    show_default=True,
    help="Limit the number of reads to classify (0 disables the limit).",
)
@click.option(
    "--max-orientation-error-1d",
    default=0.1,
    type=float,
    show_default=True,
    help="Maximum orientation error allowed for 1D classification.",
)
@click.option(
    "--max-orientation-error-2d",
    default=0.2,
    type=float,
    show_default=True,
    help="Maximum orientation error allowed for 2D classification.",
)
@click.option(
    "--min-mapped-inserts",
    default=3,
    type=int,
    show_default=True,
    help="Minimum number of mapped inserts required for classification.",
)
@click.option(
    "--backbone-prefix",
    default="BB",
    type=str,
    show_default=True,
    help="Prefix for backbone reference names to ignore during classification.",
)
@click.option(
    "--consecutiveness-method",
    default="new_occurrence",
    type=click.Choice(
        [
            "new_occurrence",
        ],
        case_sensitive=False,
    ),
    show_default=True,
    help="Method to calculate whether a segment appears consecutively in the original concatemer.",
)
def cli_consensus(
    input_bam,
    output_fastq,
    output_json,
    limit_calls,
    max_orientation_error_1d,
    max_orientation_error_2d,
    min_mapped_inserts,
    backbone_prefix,
    consecutiveness_method,
):
    click.echo("Cycas: Consensus calling on complex concatemers.")

    config = ReadClassificationConfig(
        max_orientation_error_1D=max_orientation_error_1d,
        max_orientation_error_2D=max_orientation_error_2d,
        min_mapped_inserts=min_mapped_inserts,
        backbone_prefix=backbone_prefix,
        consecutiveness_method=consecutiveness_method,
    )

    app.call_consensus(
        bam_file=input_bam,
        output_fastq=output_fastq,
        metadata_json=output_json,
        limit_calls=limit_calls,
        config=config,
    )


@cli.command(
    name="simulate",
    help="Simulate Cyclomics-seq reads from a reference genome.",
)
@click.option(
    "--reference-file",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Path to a reference genome file.",
)
@click.option(
    "--output-dir",
    required=True,
    type=click.Path(writable=True, dir_okay=True),
    help="Path to the output dir to save the results.",
)
@click.option(
    "--read-type",
    default=["1D", "3D"],
    multiple=True,
    help="Type of read to simulate.",
)
@click.option(
    "--num-reads",
    default=10_000,
    type=int,
    show_default=True,
    help="Number of reads to simulate.",
)
@click.option(
    "--insert-length",
    default=(100, 300),
    nargs=2,
    type=int,
    show_default=True,
    help="Insert length range (min, max).",
)
@click.option(
    "--num-repeats",
    default=(2, 20),
    nargs=2,
    type=int,
    show_default=True,
    help="Number of repeats (min, max).",
)
@click.option(
    "--error-rate-mismatch",
    default=0.0,
    type=float,
    show_default=True,
    help="Rate of mismatch errors.",
)
@click.option(
    "--error-rate-deletion",
    default=0.0,
    type=float,
    show_default=True,
    help="Rate of deletion errors.",
)
@click.option(
    "--error-rate-insertion",
    default=0.0,
    type=float,
    show_default=True,
    help="Rate of insertion errors.",
)
@click.option(
    "--tdt-length",
    default=(0, 0),
    nargs=2,
    type=int,
    show_default=True,
    help="Random TDT length (min, max).",
)
@click.option(
    "--num-template-switches",
    default=(0, 0),
    nargs=2,
    type=int,
    show_default=True,
    help="Number of template switches (min, max).",
)
@click.option(
    "--num-pore-chimeras",
    default=(0, 0),
    nargs=2,
    type=int,
    show_default=True,
    help="Number of pore chimeras (min, max).",
)
@click.option(
    "--num-library-chimeras",
    default=(0, 0),
    nargs=2,
    type=int,
    show_default=True,
    help="Number of library chimeras (min, max).",
)
@click.option(
    "--jagged-end-left",
    default=(0, 0),
    nargs=2,
    type=int,
    show_default=True,
    help="Jagged end length on the left (min, max).",
)
@click.option(
    "--jagged-end-right",
    default=(0, 0),
    nargs=2,
    type=int,
    show_default=True,
    help="Jagged end length on the right (min, max).",
)
@click.option(
    "--seed",
    default=0,
    type=int,
    show_default=True,
    help="Random seed for reproducibility.",
)
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    help="Overwrite output directory if it exists.",
)
def cli_simulate_reads(
    reference_file,
    output_dir,
    read_type,
    num_reads,
    insert_length,
    num_repeats,
    error_rate_mismatch,
    error_rate_deletion,
    error_rate_insertion,
    tdt_length,
    num_template_switches,
    num_pore_chimeras,
    num_library_chimeras,
    jagged_end_left,
    jagged_end_right,
    jagged_end_max_right,
    seed,
    overwrite,
):
    """
    Simulate Cyclomics-seq reads from a reference genome.
    """
    click.echo("Cycas: Cyclomics-seq read simulation.")

    if os.path.exists(output_dir):
        click.echo(f"Output directory {output_dir} already exists.")
        if overwrite:
            click.echo(f"Overwriting output directory {output_dir}.")
            shutil.rmtree(output_dir)
            os.makedirs(output_dir)
        else:
            return
    else:
        os.makedirs(output_dir)

    app.simulate_reads(
        reference_file=reference_file,
        output_dir=output_dir,
        read_type=read_type,
        num_reads=num_reads,
        insert_length=insert_length,
        num_repeats=num_repeats,
        error_rate_mismatch=error_rate_mismatch,
        error_rate_deletion=error_rate_deletion,
        error_rate_insertion=error_rate_insertion,
        tdt_length=tdt_length,
        num_template_switches=num_template_switches,
        num_pore_chimeras=num_pore_chimeras,
        num_library_chimeras=num_library_chimeras,
        jagged_end_left=jagged_end_left,
        jagged_end_right=jagged_end_right,
        seed=seed,
    )


@cli.group(name="calibration", help="Calibration subcommands")
def calibration():
    pass


@calibration.command(
    name="file-prep",
    help="Fit the calibration of a consensus calling algorithm.",
)
@click.option(
    "--hdf5-file",
    required=True,
    help="Path to the a HDF5 file with the msa data.",
    type=click.Path(exists=True),
)
@click.option(
    "--reference-file",
    required=True,
    help="Path to the reference file used in the alignment of the data",
    type=click.Path(exists=True),
)
@click.option(
    "--output-dir",
    required=True,
    help="Directory to write the output files to.",
    type=click.Path(exists=True),
)
@click.option(
    "--method",
    default="fraction",
    help="Calibration method to use.",
    type=click.Choice(
        [
            "fraction",
        ],
        case_sensitive=False,
    ),
    show_default=True,
)
@click.option(
    "--min-alignment-quality",
    default=60,
    help="Only consider alignments with this or higher alignment quality",
    show_default=True,
)
@click.option(
    "--vcf-file",
    required=False,
    help="Path to the VCF file containing variant calls.",
    type=click.Path(exists=True),
    default=None,
    show_default=True,
)
@click.option(
    "--skip-variants",
    is_flag=True,
    default=False,
    help="Skip variant calls during calibration.",
    show_default=True,
)
def cli_calibrate_consensus(
    hdf5_file,
    reference_file,
    output_dir,
    method,
    min_alignment_quality,
    vcf_file,
    skip_variants,
):
    if skip_variants and vcf_file is None:
        click.echo("VCF file must be provided with --skip-variants.")
        return None

    app.calibration_file_prep(
        hdf5_file=hdf5_file,
        reference_file=reference_file,
        output_dir=output_dir,
        method=method,
        min_alignment_quality=min_alignment_quality,
        vcf_file=vcf_file,
        skip_variants=skip_variants,
    )


@calibration.command(
    name="fit-model",
    help="Fit the model according to the prepared data.",
)
@click.option(
    "--input-dir",
    required=True,
    help="Path to the directory containing the input files.",
    type=click.Path(exists=True),
)
@click.option(
    "--output-dir",
    required=True,
    help="Directory to write the output files to.",
    type=click.Path(exists=True),
)
@click.option(
    "--method",
    default="fraction",
    help="Calibration method to use.",
    type=click.Choice(
        [
            "fraction",
        ],
        case_sensitive=False,
    ),
    show_default=True,
)
@click.option(
    "--max-depth",
    default=50,
    help="Max amount of RCA repeats to calibrate to",
    type=int,
    show_default=True,
)
@click.option(
    "--min-coverage",
    default=1000,
    help="Min amount of cases to calibrate to. For example, with 1000 we need 1000 occurrences to estimate the error of that case",
    type=int,
    show_default=True,
)
def calibrate_consensus(
    input_dir,
    output_dir,
    method,
    max_depth,
    min_coverage,
):
    app.calibration_fit_model(
        input_dir=input_dir,
        output_dir=output_dir,
        method=method,
        max_depth=max_depth,
        min_coverage=min_coverage,
    )


if __name__ == "__main__":
    cli()
