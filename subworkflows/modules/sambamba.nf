#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process SambambaSortSam{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'mgibio/dna-alignment:1.0.0'
    label 'many_med_cpu_huge_mem'

    input:
        each path(sam)

    output:
        path "*sorted.bam"
        
    script:
        """
        sambamba view -S -f bam $sam | sambamba sort -t ${task.cpus} /dev/stdin -o "${sam.simpleName}_sorted.bam"
        """
}