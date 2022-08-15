#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MinionQc{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
   
    input:
        path summary

    output:
        path "*.png" , emit: plots
        path "*.yaml" , emit: summary

    script:
        """
        $params.minionqc_location -i  $summary
        """
}

process MinionQcToJson{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        path minionqc_yaml

    output:
        tuple val("${minionqc_yaml.BaseName}"), path("global.json")

    script:
        """
        yq '.' $minionqc_yaml > global.json
        """
}
