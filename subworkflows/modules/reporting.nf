#!/usr/bin/env nextflow
nextflow.enable.dsl=2



process GenerateHtmlReport {
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container "damicyclomics/cyclomics_dash:0.2.1"
    label 'many_cpu_medium'


    input:
        tuple val(X), path(json_reads)
        tuple val(Y), path(json_globals)
        tuple val(Z), path(vcf)

    output:
        path "*.html" , emit: html

    script:
        """
        python /app/main.py \
        --sample $json_reads \
        --json_global $json_globals \
        --vcf $vcf \
        --report ${X}_${Y}_${Z}.html
        """
}


process GenerateHtmlReportWithControl {
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container "damicyclomics/cyclomics_dash:dev-unstable"
        label 'many_cpu_medium'

    input:
        tuple val(X), path(json_reads)
        tuple val(Y), path(json_globals)
        tuple val(Z), path(vcf)
        path(tvcf)

    output:
        path "*.html" , emit: html

    script:

        """
        python /app/main.py \
        --sample $json_reads \
        --json_global $json_globals \
        --vcf $vcf \
        --report ${X}_${Y}_${Z}.html \
       --vcf_truth $tvcf
        """
}
