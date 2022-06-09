#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process Cycas{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'damicyclomics/cycas:0.4.1'
    cpus = 1
    
    memory { 2.GB * task.attempt }
    time { 1.hour * task.attempt }

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    input:
        tuple val(X), path(bam), path(bai)

    output:
        tuple val(X), path("${bam.SimpleName}.consensus.fastq"), path("${bam.SimpleName}.metadata.json"), path ("plots/*")

    script:
        """
        mkdir plots
        python /app/cycas.py consensus --bam-file $bam --output ${bam.SimpleName}.consensus.fastq --metadata-json  ${bam.SimpleName}.metadata.json --plot-readtypes
        """
    }

process CycasSplit{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'damicyclomics/cycas:0.2.2-rc5'
    // label 'many_low_cpu_low_mem'
    cpus = 1
    
    memory { 2.GB * task.attempt }
    time { 1.hour * task.attempt }

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
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
