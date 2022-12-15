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
params.input_read_dir             = ""
params.read_pattern               = "{pass,fastq_pass}/**.{fq,fastq,fq.gz,fastq.gz}"
params.sequencing_quality_summary = "sequencing_summary*.txt"
params.backbone                   = "BB42"
params.backbone_name              = ""
params.region_file                = "auto"


params.reference = ""
// reference indexes are expected to be in reference folder
params.output_dir = "$HOME/Data/CyclomicsSeq"


// method selection
params.qc                   = "simple" // simple or skip
params.consensus_calling    = "cycas" // simple or skip
params.alignment            = "bwamem"  // BWA, Latal, Lastal-trained or skip
params.variant_calling      = "validate"
params.extra_haplotyping    = "skip"
params.report               = "yes"
params.split_on_adapter     = "yes"
params.quick_results        = false

// Pipeline performance metrics
params.min_repeat_count = 3

if (params.backbone == "BB41") {
    backbone_file = "$projectDir/backbones/BB41.fasta"
}
else if (params.backbone == "BB42") {
    backbone_file = "$projectDir/backbones/BB42.fasta"
}
else if (params.backbone == "BB22") {
    backbone_file = "$projectDir/backbones/BB22.fasta"
}
else if (params.backbone == "BB25") {
    backbone_file = "$projectDir/backbones/BB25.fasta"
}
else if (params.backbone == "BBCS") {
    backbone_file = "$projectDir/backbones/BBCS.fasta"
}
else if (params.backbone == "BBCR") {
    backbone_file = "$projectDir/backbones/BBCR.fasta"
}
else {
    backbone_file = params.backbone
}


// ### Printout for user
log.info """
    ===================================================
    Cyclomics/CycloSeq : Cyclomics informed pipeline
    ===================================================
    Inputs:
        input_reads              : $params.input_read_dir
        read_pattern             : $params.read_pattern
        reference                : $params.reference
        backbone                 : $params.backbone
        backbone_name            : $params.backbone_name
        region_file              : $params.region_file  
        output folder            : $params.output_dir
        Cmd line                 : $workflow.commandLine
    Method:  
        QC                       : $params.qc
        consensus_calling        : $params.consensus_calling
        Alignment                : $params.alignment
        variant_calling          : $params.variant_calling
        report                   : $params.report
        split                    : $params.split_on_adapter 

    Other:
        profile                  : $params.profile_selected
        quick_results            : $params.quick_results
"""

if (params.profile_selected == "conda"){
    log.info """
        conda environment : $params.user_conda_location
    """
}

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
    FilterWithAdapterDetection
} from "./subworkflows/filtering"

include {
    ReverseMapping
    TidehunterBackBoneQual
    CycasConsensus
    CycasMedaka
} from "./subworkflows/consensus"

include {
    Minimap2Align
    BWAAlign
    AnnotateFilter
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
    if (params.region_file != "auto"){
        if (params.variant_calling != "validate") {
            log.warn "Not All variant calling strategies support the setting of a region file!"
        }
        region_file = params.region_file
        // check if exist for fail fast behaviour
        Channel.fromPath( region_file, type: 'file', checkIfExists: true)
    }
    else {
        region_file = params.region_file
    }

    if (params.profile_selected == 'none') {
        log.warn "please set the -profile flag to `conda`, `docker` or `singularity`"
        log.warn "exiting..."
        exit(1)
    }
    if (params.profile_selected == 'local') {
        log.warn "local is available but unsupported, we advise to use a managed environment. please make sure all required software in available in the PATH"
    }
    if (params.profile_selected != 'docker') {
            if (params.consensus_calling == 'tidehunter'){
                println('Tidehunter only works with docker in this version, due to a bug in the Tidehunter release.')
                exit(1)
            }
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
    seq_summary = Channel.fromPath(sequencing_quality_summary_pattern, checkIfExists: true)
    backbone_fasta = Channel.fromPath(backbone_file, checkIfExists: true)
    
    reference_genome = Channel.fromPath(params.reference, checkIfExists: true)
    // form a pair for both .fa as well as .fasta ref genomes
    reference_genome_indexed = Channel.fromFilePairs("${params.reference}*", size: -1) { file -> file.SimpleName }
    read_info_json = ""
    PrepareGenome(reference_genome, params.reference, backbone_fasta)

    read_fastq_filtered = FilterWithAdapterDetection(read_fastq.flatten())
/*
========================================================================================
01.    Repeat identification: results in a list of read consensus in the format: val(X), path(fastq)
========================================================================================
*/
    if( params.consensus_calling == "tidehunter" ) {
        
        base_unit_reads = TidehunterBackBoneQual(read_fastq_filtered,
            reference_genome_indexed,
            backbone_fasta,
            params.tidehunter.primer_length,
            params.backbone_name,
            PrepareGenome.out.mmi_combi
        )
        read_info_json = TidehunterBackBoneQual.out.json
        base_unit_reads = TidehunterBackBoneQual.out.fastq
        split_bam = TidehunterBackBoneQual.out.split_bam
        split_bam_filtered = TidehunterBackBoneQual.out.split_bam_filtered
    }
    else if (params.consensus_calling == "cycas"){
        CycasConsensus( read_fastq_filtered,
            PrepareGenome.out.mmi_combi,
            backbone_fasta,
        )
        base_unit_reads = CycasConsensus.out.fastq
        read_info_json = CycasConsensus.out.json
        split_bam = CycasConsensus.out.split_bam
        split_bam_filtered = CycasConsensus.out.split_bam_filtered
    }
    else if(params.consensus_calling == "medaka" ) {
        CycasMedaka( read_fastq_filtered,
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
    if( params.alignment == "minimap" ) {
        Minimap2Align(base_unit_reads, PrepareGenome.out.mmi_combi, read_info_json, params.consensus_calling)
        reads_aligned = Minimap2Align.out.bam
    }
    else if( params.alignment == "bwamem" ) {
        BWAAlign(base_unit_reads, reference_genome, read_info_json, params.consensus_calling)
        reads_aligned = BWAAlign.out.bam
    }
    else if( params.alignment == "skip" ) {
        println "Skipping alignment"
    }
    else {
        error "Invalid alignment selector: ${params.alignment}"
    }
    
    // We only get the sequencing summary once we've obtained all the fastq's
    reads_aligned = AnnotateFilter(reads_aligned, seq_summary, params.min_repeat_count)

/*
========================================================================================
03.A   Variant calling
========================================================================================
*/  
    
    if( params.variant_calling == "freebayes" ) {
        FreebayesSimple(reads_aligned, PrepareGenome.out.mmi_combi)
        variant_vcf = FreebayesSimple.out
        locations = ""
    }
    else if (params.variant_calling == "varscan"){
        Varscan(reads_aligned, PrepareGenome.out.fasta_combi)
        variant_vcf = Varscan.out
        locations = ""
    }
    else if (params.variant_calling == "validate"){
        ValidatePosibleVariantLocations(
            reads_aligned,
            region_file,
            PrepareGenome.out.fasta_combi
        )
        locations = ValidatePosibleVariantLocations.out.locations
        variant_vcf = ValidatePosibleVariantLocations.out.variants
    
    }
    else {
        error "Invalid variant_calling selector: ${params.variant_calling}"
        variant_vcf = ""
        locations = ""
    }

/*
========================================================================================
04.    Reporting
========================================================================================
*/  
    PostQC(
        PrepareGenome.out.fasta_combi,
        read_fastq,
        read_fastq_filtered,
        split_bam,
        split_bam_filtered,
        base_unit_reads,
        read_info_json,
        reads_aligned,
        params.quick_results,
        locations,
        variant_vcf,
    )
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone. The results are available in following folder --> $params.output_dir\n" : "something went wrong in the pipeline" )
}
