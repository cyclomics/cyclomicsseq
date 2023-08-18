#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process Freebayes {
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'nfcore/sarek:2.7.1'
    label 'many_cpu_medium'

    input:
        tuple val(X), path(input_bam_file), path(input_bai_file), path(reference)
        tuple val(X), path(roi)

    output:
        path("${input_bam_file.SimpleName}.vcf")

    script:
        ref = reference
        """
        freebayes \
        --min-alternate-fraction 0 \
        --min-alternate-count 2 \
        --min-base-quality 10
        -t $roi \
        -f $ref \
        --vcf ${input_bam_file.SimpleName}.vcf \
        $input_bam_file
        """
}
        // ref = reference.first()