#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GenerateHtmlReport {
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container "damicyclomics/cyclomics_dash:0.0.1"

    input:
        tuple val(X), path (json)
    output:
        path "*.html" , emit: html

    script:
        """
        python /app/main.py --sample /app/tests/data/testjsons1
        """
}
