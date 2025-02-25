#!/usr/bin/env nextflow
/*
========================================================================================
    Cyclomics CyclomicsSeq pipeline
========================================================================================
    Github : https://github.com/cyclomics/cyclomicsseq-egfr
    Website: https://cyclomics.com
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    PARAMETER VALUES
========================================================================================
*/
// PARAMETERS
params.input_read_dir             = ""
params.read_pattern               = "**.{fq,fastq,fq.gz,fastq.gz}"
params.sequencing_summary_path    = "${projectDir}/sequencing_summary*.txt"
params.backbone                   = "BBCS"
params.backbone_name              = ""
params.backbone_file              = ""
params.region_file                = "auto"

params.sample_id                  = ""
params.reference                  = ""
// reference indexes are expected to be in reference folder
params.output_dir                 = "$HOME/Data/CyclomicsSeq"


// method selection
params.report                     = "detailed"
params.split_fastq_by_size        = true
params.split_on_adapter           = false
params.sequence_summary_tagging   = false
params.include_fastq_fail         = false

params.synthetics_file            = "$projectDir/contaminants/synthetic_mutants.fasta"
params.priority_limit             = (params.report == "detailed") ? 9999 : 89

// Pipeline performance metrics
params.min_repeat_count           = 3

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
    Cyclomics/CyclomicsSeq: Cyclomics informed pipeline
    ===================================================
    Inputs:
        input_read_dir           : $params.input_read_dir
        read_pattern             : $params.read_pattern
        reference                : $params.reference
        backbone                 : $params.backbone
        backbone_file            : $params.backbone_file
        region_file              : $params.region_file  
        output folder            : $params.output_dir
        Cmd line                 : $workflow.commandLine
    Method:  
        report                   : $params.report
        split                    : $params.split_on_adapter 

    Other:
        profile                  : $params.profile_selected
"""

if (params.profile_selected == "conda") {
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
    FilterWithAdapterDetection
    CycasConsensus
    AlignByID
    AlignContaminants
    PrepareGenome
    AnnotateBam
    FilterBam
    ProcessTargetRegions as ProcessContaminantRegions
    ProcessTargetRegions
    CallContaminantMutants
    CallVariantsValidate
    CrosscheckContaminantVariants
} from "./subworkflows.nf"


/*
========================================================================================
    Helper functions
========================================================================================
*/

// Function to get the first valid parent directory for a given minknow folder structure e.g.
// /data/Cyclomics/CYC000405_tRCA_test/R_n_D_sample_3/20240904_1441_P2S-02085-A_PAW62268_0eb09e0c/fastq_pass/
// 1. fastq_pass is a bad name
// 2. 20240904_1441_P2S-02085-A_PAW62268_0eb09e0c is a bad name
// 3. R_n_D_sample_3 thus is a suitable name for the sample

// if the run is barcoded and fastq_pass contains subfolders with the barcodes, those are valid sample id's!:
// ...fastq_pass/barcode03/ & ...fastq_pass/barcode04/ & ...fastq_pass/barcode05/ 
// the sample id's will be barcode03, barcode04 and barcode05

def InvalidParents = ['fastq', 'pass', 'fastq_pass', 'fastq_fail', 'fail', 'home']
// MinKnow run folder
def MinKnow_auto_run_folder_pattern = ~/^\d{8}_\d{4}_.+/  // Regular expression for 8 digits, underscore, 4 digits, underscore, and some text
def getValidParent(dir, invalidList, pattern) {
    def currentDir = dir
    while ((invalidList.contains(currentDir.simpleName) || currentDir.simpleName ==~ pattern) && currentDir.Parent != null) {
        currentDir = currentDir.Parent
    } 
    return currentDir.simpleName
}

// Used in exclusion of failed reads from the inputs.
def FastqExcludeList = ['fastq_fail', 'fail']


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
    if (params.reference == ""){
        log.warn "--reference cannot be empty!"
        exit 1
    }
    if (params.input_read_dir == ""){
        log.warn "--input_read_dir cannot be empty!"
        exit 1
    }
    // check environments
    if (params.region_file != "auto") {
        log.info "Custom region file parsing on $params.region_file"
        region_file = params.region_file
        // check if exist for fail fast behaviour
        Channel.fromPath(region_file, type: 'file', checkIfExists: true)
    }
    else {
        log.info "Auto region file parsing"
        region_file = params.region_file
    }

    if (params.profile_selected == 'local') {
        log.warn "local is available but unsupported, we advise to use a managed environment. please make sure all required software in available in the PATH"
    }

    // Process inputs:
    // add the trailing slash if its missing 
    if (params.input_read_dir.endsWith("/")) {
        read_pattern = "${params.input_read_dir}${params.read_pattern}"
    }
    else {
        read_pattern = "${params.input_read_dir}/${params.read_pattern}"
    }

    // Get the input read files, and add sample and file tags to them
    read_fastq = Channel.fromPath(read_pattern, checkIfExists: true)
    read_fastq.dump(tag: "read_fastq_all")

    // Exclude all data that is in a fastq_fail (default) or fail (rebasecalled) folder, or subfolder (e.g. barcoded runs)
    if (params.include_fastq_fail) {
        read_fastq = read_fastq
    }
    else {
        read_fastq = read_fastq.filter { 
            !(it.Parent.SimpleName in FastqExcludeList || it.Parent?.Parent.SimpleName in FastqExcludeList) 
        }
    }
    read_fastq.dump(tag: "read_fastq_no_fail")

    // Extract valid sample ID's for eah fastq file
    read_fastq = read_fastq.map { it ->
                // if parentDir is invalid, take grandParentDir etc. etc.
                def parentDir = it.Parent
                def sample_ID = getValidParent(it.Parent, InvalidParents, MinKnow_auto_run_folder_pattern)
                def file_ID = it.simpleName
                tuple(sample_ID, file_ID, it)
            }
    read_fastq.dump(tag: "read_fastq_sample_tagged")

    // Overwrite the sample ID's if provided via the cli
    if (params.sample_id != "") {
        log.warn("Overwriting sample ID for all fastq files due to cli provided sample ID")
        read_fastq = read_fastq.map { it -> tuple(params.sample_id, it[1], it[2]) }
    }
    read_fastq.dump(tag: "read_fastq")

    // Prepare reference genome, combined with backbone sequence
    backbone  = Channel.fromPath(backbone_file, checkIfExists: true)
    reference_genome = Channel.fromPath(params.reference, checkIfExists: true)
    PrepareGenome(reference_genome, params.reference, backbone)
    
    // We use .collect() to turn the genome into a value channel, enabling the input repeater
    mmi_combi = PrepareGenome.out.mmi_combi.collect()
    fasta_combi = PrepareGenome.out.fasta_combi.collect()
    bwa_index = file("${params.reference}.{,amb,ann,bwt,pac,sa}")

    seq_summary = Channel.fromPath(params.sequencing_summary_path, checkIfExists: true)

    synthetic_reads = Channel.fromPath(params.synthetics_file, checkIfExists: true)
    
    // This uses the params.split_on_adapter and params.split_fastq_by_size flags for logic control.
    read_fastq_filtered = FilterWithAdapterDetection(read_fastq)
    read_info_json = ""
/*
========================================================================================
01.    Repeat identification: results in a list of read consensus in the format: val(X), path(fastq)
========================================================================================
*/

    CycasConsensus(read_fastq_filtered, mmi_combi)
    base_unit_reads = CycasConsensus.out.fastq
    read_info_json = CycasConsensus.out.json
    split_bam = CycasConsensus.out.split_bam
    split_bam_filtered = CycasConsensus.out.split_bam_filtered

    consensus_by_id = base_unit_reads
    consensus_by_id.dump(tag: 'consensus-pre-alignment')

/*
========================================================================================
02.    Alignment
========================================================================================
*/    

    AlignByID(consensus_by_id, mmi_combi, read_info_json)
    reads_aligned = AlignByID.out

    if (params.sequence_summary_tagging) {
        reads_aligned_tagged = AnnotateBam(reads_aligned, seq_summary) 
    }
    else {
        reads_aligned_tagged = reads_aligned
    }
    reads_aligned_filtered = FilterBam(reads_aligned_tagged, params.min_repeat_count)

/*
========================================================================================
03.    Variant calling
========================================================================================
*/  
    locations = ""
    variant_vcf = ""

    ProcessTargetRegions(region_file, reads_aligned)
    regions = ProcessTargetRegions.out
    CallVariantsValidate(reads_aligned_filtered, regions, fasta_combi)
    locations = CallVariantsValidate.out.locations
    variant_vcf = CallVariantsValidate.out.variants


/*
========================================================================================
04.    Contamination check
========================================================================================
*/  
    synthetic_vcf = ""

    // String[] contaminantSeqExt = [".fasta", ".fna", ".fa", ".fastq", ".fq"]

    // if (params.synthetics_file.endsWithAny(contaminantSeqExt)) {
    //     AlignContaminants(synthetic_reads, PrepareGenome.out.fasta_ref, bwa_index)
    //     synthetics_aligned = AlignContaminants.out.bam
    //     synthetics
    //     // ProcessContaminantRegions(region_file, synthetics_aligned)
    //     ProcessContaminantRegions('auto', synthetics_aligned)

    //     CallContaminantMutants(synthetics_aligned,  ProcessContaminantRegions.out, PrepareGenome.out.fasta_combi)
    //     synthetic_vcf = CallContaminantMutants.out.variants       

    // } else if (params.synthetics_file.endsWith(".vcf")) {
    //     synthetic_vcf = synthetic_reads

    // } else {
    //     error "Invalid --synthetics_file format: ${params.synthetics_file}"
    // }

    if (synthetic_vcf) {
        contaminants_vcf = CrosscheckContaminantVariants(variant_vcf, synthetic_vcf)
    }
    else {
        contaminants_vcf = ""
    }
/*
========================================================================================
05.    Reporting
========================================================================================
*/  
    if (params.report in ["detailed", "standard"]) {
        Report(
            fasta_combi,
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
            contaminants_vcf,
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
	log.info(workflow.success ? "\nDone. The results are available in following folder --> $params.output_dir\n" : "something went wrong in the pipeline")
}
