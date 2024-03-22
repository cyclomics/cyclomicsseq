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
params.output_dir    = "CyclomicsSeq_${workflow.revision}_${workflow.runName}_${workflow.sessionId}"


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


include {
    FilterShortReads
} from "./nextflow_utils/sequence_analysis/modules/seqkit"

workflow {
     log.info """Cyclomics consensus pipeline started"""

    // Process inputs:
    // add the trailing slash if its missing 
    if (params.input_read_dir.endsWith("/")){
        read_pattern = "${params.input_read_dir}${params.read_pattern}"
    }
    else {
        read_pattern = "${params.input_read_dir}/${params.read_pattern}"
    }
    read_dir_ch = Channel.fromPath( params.input_read_dir, type: 'dir', checkIfExists: true)
    // Create an item where we have the path and the sample ID.
    //  If the path is a directory with fastq's in it directly,
    // The sample ID will be that directory name. eg. fastq_pass/ could be the sample ID.
    // Alternatively, the sample ID could be barcode01/barcode02 etc.
    reads = Channel.fromPath(read_pattern, checkIfExists: true) \
        | map(x -> [x.Parent.simpleName, x.simpleName,x])
    if (params.devtooling.show_data_structures) {reads.view()}
    
    reads_filtered = FilterShortReads(reads)
    if (params.devtooling.show_data_structures) {reads_filtered.view()}

    consensus = CygnusConsensusGeneration()
}