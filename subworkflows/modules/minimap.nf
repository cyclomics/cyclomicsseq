#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process MinimapAlign{
    // Sams are large and should never be the final data point
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'staphb/minimap2:2.23'
    label 'few_memory_intensive'

    input:
        tuple val(X), path(fasta), path(reference_genome)

    output:
        tuple val(X), path("${X}.sam") 

    script:
        """
        minimap2 -ax map-ont -t ${task.cpus} $reference_genome $fasta > ${X}.sam
        """
}

process MinimapAlignMany{
    // Same 
    // Sams are large and should never be the final data point
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'staphb/minimap2:2.23'
    label 'minimap_large'

    input:
        each path(fasta)
        path(reference_genome)

    output:
        tuple val("${fasta.baseName}"), path("${fasta.baseName}.sam") 

    script:
        """
        minimap2 -ax map-ont -t ${task.cpus} $reference_genome $fasta > ${fasta.baseName}.sam
        """
}

process Minimap2Index{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'staphb/minimap2:2.23'
    label 'few_memory_intensive'

     input:
        path(reference_genome)
    
    output:
        path("${reference_genome.baseName}.mmi")

    script:
        """
        minimap2  -ax map-ont -t ${task.cpus} -d ${reference_genome.baseName}.mmi $reference_genome
        """
}