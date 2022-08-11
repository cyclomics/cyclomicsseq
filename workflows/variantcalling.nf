#!/usr/bin/env nextflow
/*
========================================================================================
    Cyclomics/VariantCalling
========================================================================================
    Github : tbd
    Website: https://cyclomics.com
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

params.input_bam_files = "/media/dami/cyclomics_003/raw_data/ENZA_0003/consensus/align/B-8236-001.rotated.bam"
params.reference = "/home/dami/projects/ENZA005_real_data_analysis/reference/Viroflay_ref_F7-R8.fasta"
params.output_dir = "data/$workflow.runName"

// ### Printout for user
log.info """
    ===================================================
    Cyclomics/variantcalling: call freebayes on bams
    ===================================================
        input_bam_files : $params.input_bam_files
        reference       :$params.reference
        output_dir      :$params.output_dir
"""

process SamtoolsFaidx{
    // Given a sam or bam file, make a sorted bam
    // Does break when sam is not `propper` eg no @SQ tags
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'biocontainers/samtools:v1.7.0_cv4'

    cpus 1

    input:
        path(fasta)

    output:
        path "${fasta}.fai"

    script:
        """
        samtools faidx $fasta
        """
}

process Freebayes {
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'nfcore/sarek:2.7.1'
    cpus 1
    input:
        each path(input_bam_file)
        path(reference)

    output:
        path("${input_bam_file.SimpleName}.vcf")

    script:
        ref = reference.first()
        """
        freebayes -f $ref $input_bam_file -v ${input_bam_file.SimpleName}.vcf
        """
}




workflow{
    bams = Channel.fromPath(params.input_bam_files, checkIfExists: true)
    reference = Channel.fromPath(params.reference, checkIfExists: true)
    faidx = SamtoolsFaidx(reference)
    Freebayes(bams, reference.combine(faidx))

}