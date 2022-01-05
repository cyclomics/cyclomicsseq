
nextflow.enable.dsl=2

process MinionQc{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'evolbioinfo/minionqc:v1.4.1'

    input:
        path summary

    output:
        path "*.html" , emit: html
        path "*.json" , emit: json

    script:
        """
        Rscript MinIONQC.R -i  $summary/sequencing_summary.txt
        """
}

process MinionQcToJson {
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'stedolan/jq'

    input:
        path minionqc_yaml

    output:
        tuple val(X), path("QC.json")

    script:
        """
        """
}