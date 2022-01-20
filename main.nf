#!/usr/bin/env nextflow
/*
========================================================================================
    Cyclomics/CycloSeq Informed pipeline
========================================================================================
    Github : https://github.com/cyclomics/cycloseq
    Website: https://cyclomics.com
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    PARAMETER VALUES
========================================================================================
*/
// ### PARAMETERS
params.input_read_dir             = "$HOME/data/raw_data/MAR6252/"
params.read_pattern               = "fastq_pass/**.{fq,fastq}"
params.sequencing_quality_summary = "sequencing_summary*.txt"
params.backbone_fasta             = "$HOME/data/backbones/backbones/backbones_db_valid.fasta"
params.backbone_name              = "BB22"


params.reference = "$HOME/data/references/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa"
// reference indexes are expected to be in reference folder
params.output_dir = "$HOME/data/nextflow/Cyclomics_informed"

// ### Printout for user
log.info """
    ===================================================
    Cyclomics/CycloSeq : Cyclomics informed pipeline
    ===================================================
        input_reads   : $params.input_read_dir
        read_pattern  : $params.read_pattern
        reference     : $params.reference
        backbone_fasta: $params.backbone_fasta
        backbone_name : $params.backbone_name

        output folder : $params.output_dir
"""

/*
========================================================================================
    Include statements
========================================================================================
*/
include {
    QC_MinionQc
} from "./subworkflows/QC"

include {
    ReverseMapping
    TidehunterBackBoneQual
    ConsensusTroughAlignment
} from "./subworkflows/consensus"

include {
    Minimap2Align
} from "./subworkflows/align"

include {
    FreebayesSimple
    Mutect2
} from "./subworkflows/variant_calling"

include {
    Report
} from "./subworkflows/report"

/*
========================================================================================
    Workflow
========================================================================================
*/
workflow {
    read_pattern = "${params.input_read_dir}${params.read_pattern}"
    sequencing_quality_summary_pattern = "${params.input_read_dir}/${params.sequencing_quality_summary}"

    read_dir_ch = Channel.fromPath( params.input_read_dir, type: 'dir', checkIfExists: true)
    read_fastq = Channel.fromPath(read_pattern, checkIfExists: true)
    // read_fastq.view()
    qc_seq_file_ch = Channel.fromPath(sequencing_quality_summary_pattern, checkIfExists: false)
    backbone_fasta = Channel.fromPath(params.backbone_fasta, checkIfExists: true)
    
    reference_genome_raw = Channel.fromPath(params.reference, checkIfExists: true)
    // form a pair for both .fa as well as .fasta ref genomes
    reference_genome_indexed = Channel.fromFilePairs("${params.reference}*", size: -1) { file -> file.SimpleName }


    QC_MinionQc(params.input_read_dir)
/*
========================================================================================
01.    Repeat identification: results in a list of read consensus in the format: val(X), path(fastq)
========================================================================================
*/

    // ReverseMapping(read_fastq,backbone_fasta)
    TidehunterBackBoneQual(read_fastq.flatten(),
        reference_genome_indexed,
        backbone_fasta,
        params.tidehunter.primer_length,
        params.backbone_name
    )

/*
========================================================================================
02.    Alignment
========================================================================================
*/
    // TODO: Better naming of the merged bam (eg remove number extension)
    // TODO: Annotate the info from the sequencing summary eg: AnnotateSequencingSummary(Minimap2Align.out, )
    // TODO: Annotate the info from the Tidehunter summary eg: AnnotateTidehunterSummary(Minimap2Align.out, )

    Minimap2Align(TidehunterBackBoneQual.out.fastq, reference_genome_raw)
/*
========================================================================================
03.    Variant calling
========================================================================================
*/
    FreebayesSimple(Minimap2Align.out, reference_genome_raw)
    // Mutect2(Minimap2Align.out, reference_genome_raw)

/*
========================================================================================
03.    Reporting
========================================================================================
*/  

    Report(TidehunterBackBoneQual.out.json, 
        QC_MinionQc.out, 
        FreebayesSimple.out
    )
}
