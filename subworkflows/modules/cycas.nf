#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process Cycas{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/consensus", mode: 'copy'

    label 'many_low_cpu_tiny_mem'

    input:
        tuple val(X), path(bam), path(bai)

    output:
        tuple val(X), path("${bam.SimpleName}.consensus.fastq"), path("${bam.SimpleName}.metadata.json")

    script:
        
        """
        mkdir plots
        python $params.cycas_location consensus --bam-file $bam --output ${bam.SimpleName}.consensus.fastq --metadata-json  ${bam.SimpleName}.metadata.json
        """
    }

process CycasSplit{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'damicyclomics/cycas:0.2.2-rc5'
    label 'many_low_cpu_low_mem'

    input:
        tuple val(X), path(bam), path(bai)

    output:
        path("${X}.tar")

    script:
        """
        mkdir /cycas_results
        python /app/cycas.py split --bam-file $bam --output /cycas_results/${X} --create-metadata-json --create-classification-detail-json
        tar -cvf ${X}.tar /cycas_results 
        """
    }
