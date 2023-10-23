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
params.read_pattern               = "**.{fq,fastq,fq.gz,fastq.gz}"
params.sequencing_summary_path = "${projectDir}/sequencing_summary*.txt"
params.backbone                   = "BBCS"
params.backbone_name              = ""
params.region_file                = "auto"


params.reference = ""
// reference indexes are expected to be in reference folder
params.output_dir = "$HOME/Data/CyclomicsSeq"


// method selection
params.report                   = "detailed"
params.consensus_calling        = "cycas"
params.alignment                = "bwamem"
params.variant_calling          = "validate"
params.report                   = true
params.split_on_adapter         = false
params.sequence_summary_tagging = false
params.backbone_file            = ""
params.priority_limit = (params.report == "detailed") ? 9999 : 89

// Pipeline performance metrics
params.min_repeat_count = 3

if (params.backbone == "BB41") {
    backbone_file = "$projectDir/backbones/BB41.fasta"
}
else if (params.backbone == "BB41T") {
    backbone_file = "$projectDir/backbones/BB41T.fasta"
}
else if (params.backbone == "BB42") {
    backbone_file = "$projectDir/backbones/BB42.fasta"
}
else if (params.backbone == "BB43") {
    backbone_file = "$projectDir/backbones/BB43.fasta"
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
    backbone_file = params.backbone_file
}


// ### Printout for user
log.info """
    ===================================================
    Cyclomics/CyclomicsSeq : Cyclomics informed pipeline
    ===================================================
    Inputs:
        input_reads              : $params.input_read_dir
        read_pattern             : $params.read_pattern
        reference                : $params.reference
        backbone                 : $params.backbone
        backbone_file            : $params.backbone_file
        region_file              : $params.region_file  
        output folder            : $params.output_dir
        Cmd line                 : $workflow.commandLine
    Method:  
        report                   : $params.report
        consensus_calling        : $params.consensus_calling
        Alignment                : $params.alignment
        variant_calling          : $params.variant_calling
        report                   : $params.report
        split                    : $params.split_on_adapter 

    Other:
        profile                  : $params.profile_selected
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
    Report
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
    AnnotateBam
    FilterBam
    PrepareGenome
} from "./subworkflows/align"

include {
    ProcessTargetRegions
    CallVariantsFreebayes
    Mutect2
    ValidatePosibleVariantLocations
} from "./subworkflows/variant_calling"

// include {
//     Report
// } from "./subworkflows/report"

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

    // if (params.profile_selected == 'none') {
    //     log.warn "please set the -profile flag to `conda`, `docker` or `singularity`"
    //     log.warn "exiting..."
    //     exit(1)
    // }
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

    read_dir_ch = Channel.fromPath( params.input_read_dir, type: 'dir', checkIfExists: true)
    read_fastq = Channel.fromPath(read_pattern, checkIfExists: true)
    seq_summary = Channel.fromPath(params.sequencing_summary_path, checkIfExists: true)
    backbone_fasta = Channel.fromPath(backbone_file, checkIfExists: true)
    
    reference_genome = Channel.fromPath(params.reference, checkIfExists: true)
    // form a pair for both .fa as well as .fasta ref genomes
    reference_genome_indexed = Channel.fromFilePairs("${params.reference}*", size: -1) { file -> file.SimpleName }
    bwa_index = file( "${params.reference}.{,amb,ann,bwt,pac,sa}")

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
        BWAAlign(base_unit_reads, PrepareGenome.out.fasta_ref , bwa_index, read_info_json)
        reads_aligned = BWAAlign.out.bam
    }
    else if( params.alignment == "skip" ) {
        println "Skipping alignment"
    }
    else {
        error "Invalid alignment selector: ${params.alignment}"
    }
    
    // We only get the sequencing summary once we've obtained all the fastq's
    if (params.sequence_summary_tagging) {
        reads_aligned_tagged = AnnotateBam(reads_aligned, seq_summary) 
    }
    else {
        reads_aligned_tagged = reads_aligned
    }
    reads_aligned_filtered = FilterBam(reads_aligned_tagged, params.min_repeat_count)
/*
========================================================================================
03.A   Variant calling
========================================================================================
*/  
    ProcessTargetRegions(region_file, reads_aligned_filtered)
    regions = ProcessTargetRegions.out
    locations = ""
    variant_vcf = ""

    if (params.variant_calling == "validate") {
        ValidatePosibleVariantLocations(
            reads_aligned_filtered,
            regions,
            PrepareGenome.out.fasta_combi
        )
        locations = ValidatePosibleVariantLocations.out.locations
        variant_vcf = ValidatePosibleVariantLocations.out.variants
    }
    else if (params.variant_calling == "freebayes") {
        CallVariantsFreebayes(
            reads_aligned_filtered,
            regions,
            PrepareGenome.out.fasta_combi
        )
        locations = CallVariantsFreebayes.out.locations
        variant_vcf = CallVariantsFreebayes.out.variants
    }
    else {
        error "Invalid variant_calling selector: ${params.variant_calling}"
    }

/*
========================================================================================
04.    Reporting
========================================================================================
*/  
    if (params.report in ["detailed", "standard"]){
        Report(
            PrepareGenome.out.fasta_combi,
            read_fastq,
            read_fastq_filtered,
            split_bam,
            split_bam_filtered,
            base_unit_reads,
            read_info_json,
            reads_aligned_filtered,
            regions,
            locations,
            variant_vcf,
        )
    }
    else if (params.report == "skip") {
        println('Skipping report generation entirely')
    }
    else {
        error "Invalid reporting selector: ${params.report}"
    }

}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone. The results are available in following folder --> $params.output_dir\n" : "something went wrong in the pipeline" )
}
