#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process PerbaseBaseDepth {
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'few_very_memory_intensive'

    input:
        tuple val(X), path(input_bam_file), path(input_bai_file), path(reference)
        tuple val(X_bed), path(bed)
        val(output_name)

    output:
        tuple val(X), path(output_name)

    script:
    """
    samtools faidx $reference
    perbase base-depth -F 256 -D $params.perbase.max_depth -t $task.cpus -z -b $bed --ref-fasta $reference $input_bam_file > $output_name
    """
}
