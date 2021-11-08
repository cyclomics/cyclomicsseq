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
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/




//  main channels with the input dir and the reads
// read_dir_ch = Channel.fromPath( params.input_read_dir, type: 'dir', checkIfExists: true)
// read_file_ch = Channel.fromPath(read_pattern, checkIfExists: true)

// // qc file
// qc_seq_file_ch = Channel.create()
// Channel.fromPath(sequencing_quality_summary_pattern).into(qc_seq_file_ch)
include {BwaIndex;
        BwaMemSorted
} from "./modules/bwa.nf"

include {
    QC_pycoqc;
} from "./subworkflows/QC.nf"

include {
    FilteringBasic;
} from "./subworkflows/filtering.nf"

include {
    TideHunterBasic;
} from "./subworkflows/repeat_identifier.nf"


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
    01.    Quality Control
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
    01.    Filtering
    ========================================================================================
    */ 
    if( params.filtering == 'simple' )
        repeats = FilteringBasic(read_file_ch.flatten())
    else if( params.qc == 'skip' )
        println "Skipping  Filtering"
    else
        error "Invalid repeat finder selector: ${params.filtering}"

    // /*
    // ========================================================================================
    // 02.    repeat splitting
    // ========================================================================================
    // */ 
    if( params.repeat_split == 'simple' )
        repeats = TideHunterBasic(read_file_ch)
    else if( params.repeat_split == 'skip' )
        println "Skipping  Filtering"
    else
        error "Invalid repeat finder selector: ${params.repeat_split}"


    reference_tmp = Channel.fromFilePairs("${params.reference}*.{fa,fasta}*", size: -1).first()

    println "${read_file_ch.view()}"
    println "${reference_tmp.view()}"
    
    BwaMemSorted(
        read_file_ch,
        reference_genome_indexed,
    )

    // if( params.repeat_splitting == 'simple' )
    //     repeats = RepeatSplitBasic(input_reads, input_reads)
    // else
    //     error "Invalid repeat finder selector: ${params.repeat_splitting}"
    // /*
    // ========================================================================================
    // 04.    Consensus calling Identification
    // ========================================================================================
    // */ 
    // if( params.consensus_calling == 'simple' )
    //     repeats = ConsensusBasic(input_reads, input_reads)
    // else
    //     error "Invalid repeat finder selector: ${params.consensus_calling}"
    }



























    // Extract3PrimeFasta(params.backbone_fasta, params.backbone_prime_length)
    // Extract5PrimeFasta(params.backbone_fasta, params.backbone_prime_length)

    // Tidehunter(input_reads)
    // TidehunterFullLength(input_reads,
    //     Extract3PrimeFasta.out,
    //     Extract5PrimeFasta.out)

    // TidehunterAggregate(Tidehunter.out,
    // TidehunterFullLength.out)
    // TideHunterTrimmmerPrimer(TidehunterAggregate.out)

    // Tidehunter53(
    //     input_reads,
    //     Extract3PrimeFasta.out,
    //     Extract5PrimeFasta.out
    // )
//     Tidehunter(
//         input_reads
//     )

//     // BwaIndex(
//     //     input_reads_count
//     // )
// // align consensus to the expected target
//     BwaMemSorted(
//         Tidehunter.out.consensus,
//         indexed_ref,
//     )
// }

/*
========================================================================================
    THE END
========================================================================================
*/
