#!/usr/bin/env nextflow
/*
========================================================================================
    Cyclomics/CycloSeq
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
params.input_read_dir             = "/home/dami/data/raw_data/MAR6228/"
params.read_pattern               = "fastq_pass/**.{fq,fastq}"
params.sequencing_quality_summary = "sequencing_summary*.txt"

params.reference = "/home/dami/data/references/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/"
params.reference_indexes = 

params.qc                   = "simple"
params.filtering            = "simple"
params.repeat_split         = "simple"
params.repeat_splitting     = "simple"
params.consensus_calling    = "simple"

/*
========================================================================================
    Reference subworkflows
========================================================================================
*/

include {
    QC_pycoqc;
} from "./subworkflows/QC.nf"

include {
    FilteringBasic;
} from "./subworkflows/filtering.nf"

include {
    TideHunterBasic;
} from "./subworkflows/repeat_identifier.nf"

include {
    AlignBWA as Align
} from "./subworkflows/align.nf"

workflow {
    read_pattern = "${params.input_read_dir}${params.read_pattern}"
    sequencing_quality_summary_pattern = "${params.input_read_dir}/${params.sequencing_quality_summary}"

    read_dir_ch = Channel.fromPath( params.input_read_dir, type: 'dir', checkIfExists: true)
    read_file_ch = Channel.fromPath(read_pattern, checkIfExists: true)

    qc_seq_file_ch = Channel.fromPath(sequencing_quality_summary_pattern, checkIfExists: false)
    
    // form a pair for both .fa as well as .fasta ref genomes
    reference_genome_indexed = Channel.fromFilePairs("${params.reference}*.{fa,fasta}*", size: -1)
    // TODO: check for .amb,.ann,.bwt,.pac,.sa
    // TODO: index when no index files present


    println "Starting workflow"

    /*
    ========================================================================================
    00.    Quality Control
    ========================================================================================
    */
    if( params.qc == 'simple' )
        qc_report = QC_pycoqc(qc_seq_file_ch)

    else if( params.qc == 'skip' )
        println "Skipping QC control"

    else
        error "Invalid qc selector: ${params.qc}"

    /*
    ========================================================================================
    01.    repeat splitting
    ========================================================================================
    */ 
    if( params.repeat_split == 'simple' )
        repeats = TideHunterBasic(read_file_ch, reference_genome_indexed)
    else if( params.repeat_split == 'skip' )
        println "Skipping  Filtering"
    else
        error "Invalid repeat finder selector: ${params.repeat_split}"

    /*
    ========================================================================================
    02.    Filtering
    ========================================================================================
    */ 
    if( params.filtering == 'simple' )
        println "Not implemented"
        exit 1
        filtered_reads = FilteringBasic(read_file_ch.flatten())
    else if( params.qc == 'skip' )
        filtered_reads = read_file_ch
        println "Skipping  Filtering"
    else
        error "Invalid repeat finder selector: ${params.filtering}"

     /*
    ========================================================================================
    03.    Align
    ========================================================================================
    */ 
    Align(filtered_reads)