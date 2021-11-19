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

params.reference = "/home/dami/data/references/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa"
// freference indexes are expected to be in reference folder
params.output_dir = "data/out/$workflow.runName"

params.qc                   = "simple" // simple or skip
params.filtering            = "simple" // simple or skip
params.repeat_splitting     = "simple" // simple, longest, tidehunterquality or skip
params.consensus_calling    = "simple" // simple or skip
params.final_alignment      = "BWA"  // BWA, Latal, Lastal-trained or skip

// ### Printout for user
log.info """
    ===================================================
    Cyclomics/CycloSeq : Cyclomics consensus pipeline
    ===================================================
        input_reads   : $params.input_read_dir
        read_pattern  : $params.read_pattern
        reference     : $params.reference
"""

/*
========================================================================================
    Loading
========================================================================================
*/

include {
    QC_pycoqc;
} from "./subworkflows/QC"

include {
    FilteringBasic
} from "./subworkflows/filtering"

include {
    TideHunterBasic
    TideHunterKeepLongest
} from "./subworkflows/repeat_identifier"

include {
    AlignBWA
} from "./subworkflows/align"

include {
    PostAlignmentQC
} from "./subworkflows/post_align_qc"

workflow {
    read_pattern = "${params.input_read_dir}${params.read_pattern}"
    sequencing_quality_summary_pattern = "${params.input_read_dir}/${params.sequencing_quality_summary}"

    read_dir_ch = Channel.fromPath( params.input_read_dir, type: 'dir', checkIfExists: true)
    read_file_ch = Channel.fromPath(read_pattern, checkIfExists: true)

    qc_seq_file_ch = Channel.fromPath(sequencing_quality_summary_pattern, checkIfExists: false)
    
    // form a pair for both .fa as well as .fasta ref genomes
    reference_genome_indexed = Channel.fromFilePairs("${params.reference}*", size: -1) { file -> file.SimpleName }
    // println reference_genome_indexed.view()
    // TODO: check for .amb,.ann,.bwt,.pac,.sa
    // TODO: index when no index files present


    println "Starting workflow"

    /*
    ========================================================================================
    00.    Quality Control
    ========================================================================================
    */

    if( params.qc == "simple" )
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

    if( params.repeat_splitting == "simple" )
        repeats = TideHunterBasic(read_file_ch, reference_genome_indexed)
    
    else if( params.repeat_splitting == 'longest')
        repeats = TideHunterKeepLongest(read_file_ch, reference_genome_indexed)
    
    else if( params.repeat_splitting == 'tidehunterquality')
        repeats = TideHunterQuality(read_file_ch, reference_genome_indexed)

    else if( params.repeat_splitting == 'skip' )
        println "Skipping  Filtering"

    else
        error "Invalid repeat finder selector: ${params.repeat_splitting}"

    /*
    ========================================================================================
    02.    Filtering
    ========================================================================================
    */ 
    // if( params.filtering == "simple" )
    //     reads_filtered = FilteringBasic(repeats.flatten())
    // else if( params.filtering == 'skip' ) {
    //     reads_filtered = repeats
    //     println "Skipping  Filtering"
    // }
    // else
    //     error "Invalid repeat finder selector: ${params.filtering}"
    


    /*
    ========================================================================================
    03.    Alignment
    ========================================================================================
    */ 

    // dont call .out on assigned workflow outputs
    // AlignBWA(repeats,  reference_genome_indexed)

    PostAlignmentQC(repeats)
}
