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
params.input_read_dir             = "$HOME/Data/raw_data/MAR6252/"
params.read_pattern               = "fastq_pass/**.{fq,fastq,fq.gz,fastq.gz}"
params.sequencing_quality_summary = "sequencing_summary*.txt"
params.backbone_fasta             = "$HOME/Data/backbones/backbones_db_current_slim.fasta"
params.backbone_name              = ""

params.validation_location_file   = ""

params.reference = "$HOME/Data/references/Homo_sapiens/T2T/chm13v2.0.fa"
// reference indexes are expected to be in reference folder
params.output_dir = "$HOME/Data/CyclomicsSeq"


// method selection
params.qc                   = "simple" // simple or skip
params.consensus_calling    = "cycas" // simple or skip
params.alignment            = "minimap"  // BWA, Latal, Lastal-trained or skip
params.variant_calling      = "varscan"
params.extra_haplotyping    = "skip"
params.report               = "yes"
params.quick_results        = false



// ### Printout for user
log.info """
    ===================================================
    Cyclomics/CycloSeq : Cyclomics informed pipeline
    ===================================================
    Inputs:
        input_reads              : $params.input_read_dir
        read_pattern             : $params.read_pattern
        reference                : $params.reference
        backbone_fasta           : $params.backbone_fasta
        backbone_name            : $params.backbone_name
        validation_location_file : $params.validation_location_file  
        output folder            : $params.output_dir
    Method:  
        QC                       : $params.qc
        consensus_calling        : $params.consensus_calling
        Alignment                : $params.alignment
        variant_calling          : $params.variant_calling
        report                   : $params.report

    Other:
        profile                  : $params.profile_selected
        quick_results            : $params.quick_results
"""

/*
========================================================================================
    Include statements
========================================================================================
*/
include {
    QC_MinionQc
    PostQC
} from "./subworkflows/QC"

include {
    ReverseMapping
    TidehunterBackBoneQual
    CycasConsensus
    CycasMedaka
} from "./subworkflows/consensus"

include {
    Minimap2Align
    Annotate
    PrepareGenome
} from "./subworkflows/align"

include {
    FreebayesSimple
    Mutect2
    Varscan
    ValidatePosibleVariantLocations
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
    /*
    ========================================================================================
    AA. Parameter processing 
    ========================================================================================
    */
    // check environments
    if (params.validation_location_file != ""){
        log.warn "Overwriting variant_calling strategy due to the presence of a --validation_location_file"
        params.variant_calling = "validate"
        log.info "--variant_calling set to $params.variant_calling"
    }
    if (params.profile_selected == 'none') {
        log.warn "please set the -profile flag to `conda`, `docker` or `singularity`"
        log.warn "exiting..."
        exit(1)
    }
    if (params.profile_selected == 'local') {
        log.warn "local is available but unsupported, we advise to use a managed environment. please make sure all required software in available in the PATH"
    }

    // Process inputs:
    // add the trailing slash if its missing 
    if (params.input_read_dir.endsWith("/")){
        read_pattern = "${params.input_read_dir}${params.read_pattern}"
    }
    else {
        read_pattern = "${params.input_read_dir}/${params.read_pattern}"
    }

    sequencing_quality_summary_pattern = "${params.input_read_dir}/${params.sequencing_quality_summary}"

    read_dir_ch = Channel.fromPath( params.input_read_dir, type: 'dir', checkIfExists: true)
    read_fastq = Channel.fromPath(read_pattern, checkIfExists: true)
    // read_fastq.view()
    seq_summary = Channel.fromPath(sequencing_quality_summary_pattern, checkIfExists: false)
    backbone_fasta = Channel.fromPath(params.backbone_fasta, checkIfExists: true)
    
    reference_genome = Channel.fromPath(params.reference, checkIfExists: true)
    // form a pair for both .fa as well as .fasta ref genomes
    reference_genome_indexed = Channel.fromFilePairs("${params.reference}*", size: -1) { file -> file.SimpleName }
    read_info_json = ""
    
    PrepareGenome(reference_genome, params.reference, backbone_fasta)

/*
========================================================================================
00.    raw data quality control
========================================================================================
*/
    if( params.qc == "simple" ) {

        QC_MinionQc(seq_summary)
    }
    else if( params.qc == "skip" ) {
        println "Skipping QC control"
    }
    else {
        error "Invalid qc selector: ${params.qc}"
    }

/*
========================================================================================
01.    Repeat identification: results in a list of read consensus in the format: val(X), path(fastq)
========================================================================================
*/
    // ReverseMapping(read_fastq,backbone_fasta)h
    if( params.consensus_calling == "tidehunter" ) {
        base_unit_reads = TidehunterBackBoneQual(read_fastq.flatten(),
            reference_genome_indexed,
            backbone_fasta,
            params.tidehunter.primer_length,
            params.backbone_name
        )
        read_info_json = TidehunterBackBoneQual.out.json
        base_unit_reads = TidehunterBackBoneQual.out.fastq
    }
    else if (params.consensus_calling == "cycas"){
        CycasConsensus( read_fastq.flatten(),
            PrepareGenome.out.mmi_combi,
            backbone_fasta,
        )
        base_unit_reads = CycasConsensus.out.fastq
        read_info_json = CycasConsensus.out.json
    }
    else if(params.consensus_calling == "medaka" ) {
        CycasMedaka( read_fastq.flatten(),
            PrepareGenome.out.mmi_combi,
            backbone_fasta,
        )
        base_unit_reads = CycasMedaka.out.fastq
        read_info_json = CycasMedaka.out.json
        // reference_genome_raw = 
    }
    else if(params.consensus_calling == "skip" ) {
        println "Skipping consensus_calling"
    }
    else {
        error "Invalid consensus_calling selector: ${params.consensus_calling}"
    }
/*
========================================================================================
02.    Alignment
========================================================================================
*/
    // TODO: Annotate the info from the Tidehunter summary eg: AnnotateTidehunterSummary(Minimap2Align.out, )
    
    if( params.alignment == "minimap" ) {
        Minimap2Align(base_unit_reads, PrepareGenome.out.mmi_combi)
        reads_aligned = Minimap2Align.out.bam
    }
    else if( params.alignment == "skip" ) {
        println "Skipping alignment"
    }
    else {
        error "Invalid alignment selector: ${params.alignment}"
    }
    
    // We only get the sequencing summary once we've obtained all the fastq's
    reads_aligned = Annotate(reads_aligned, seq_summary)

/*
========================================================================================
03.A   Variant calling
========================================================================================
*/  
    
    if( params.variant_calling == "freebayes" ) {
        FreebayesSimple(reads_aligned, PrepareGenome.out.mmi_combi)
        vcf = FreebayesSimple.out
    }
    else if (params.variant_calling == "varscan"){
        Varscan(reads_aligned, PrepareGenome.out.fasta_combi)
        vcf = Varscan.out
    }
    else if (params.variant_calling == "validate"){
        ValidatePosibleVariantLocations(
            reads_aligned,
            params.validation_location_file,
            PrepareGenome.out.fasta_combi
        )
    }
    else {
        error "Invalid variant_calling selector: ${params.variant_calling}"
        vcf = ""
    }
    // else if( params.variant_calling == "skip" ) {
    //     println "Skipping variant_calling"
    //     vcf = ""
    // }
    // 

/*
========================================================================================
04.    Reporting
========================================================================================
*/ 
    PostQC(
        read_info_json,
        read_fastq,
        base_unit_reads,
        reads_aligned,
        params.quick_results,
        CycasConsensus.out.split_bam
    )
}
/*
========================================================================================
04.    QC stuff
========================================================================================
*/  

    // Count reads in input fastqs

    // Count classifications
