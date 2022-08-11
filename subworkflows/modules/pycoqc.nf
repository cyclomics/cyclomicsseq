#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process PYCOQC {
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container "quay.io/biocontainers/pycoqc:2.5.2--py_0"
    label 'many_cpu_medium'


    input:
        path summary

    output:
        path "*.html" , emit: html
        path "*.json" , emit: json

    script:
    """
    pycoQC \\
        -f $summary \\
        -o pycoqc.html \\
        -j pycoqc.json
    """
}
