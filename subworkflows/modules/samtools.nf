#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process SamtoolsIndex{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'biocontainers/samtools:v1.7.0_cv4'
    label 'many_cpu_medium'


    input:
        each path(bam)

    output:
        tuple path(bam), path("*.bai")

    script:
        """
        samtools index $bam
        """
}

process SamtoolsSort{
    // Given a sam or bam file, make a sorted bam
    // Does break when sam is not `propper` eg no @SQ tags
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'biocontainers/samtools:v1.7.0_cv4'
    label 'many_cpu_medium'

    input:
        path(input_sam)

    output:
        path "${input_sam.SimpleName}_sorted.bam"

    script:
        """
        samtools sort $input_sam -O bam -o "${input_sam.SimpleName}_sorted.bam"
        """
}

process SamtoolsDepth{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'biocontainers/samtools:v1.7.0_cv4'
    label 'many_cpu_medium'

    input:
        tuple val(X), path(input_bam), path(bai)

    output:
        path "${input_bam.SimpleName}.bed"

    script:
        """
        samtools depth $input_bam > ${input_bam.SimpleName}.bed
        """
}

process SamtoolsDepthToJson{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    // the slim version of buster is missing ps, which is needed for nextflow
    container 'python:3.9.10-buster'
    label 'many_cpu_medium'

    input:
        path(input_bed)

    output:
        path "${input_bed.SimpleName}.json"

    script:
        """
        bed_to_json.py --bed $input_bed
        """
}


process SamToBam{
    // Sort, convert and index 
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'biocontainers/samtools:v1.7.0_cv4'
    label 'many_cpu_medium'

    input:
        tuple val(X), path(sam_file)

    output:
        tuple val(X), path("${X}.bam"), path("${X}.bam.bai")

    script:
        """
        samtools view -Sb $sam_file | samtools sort /dev/stdin -o ${X}.bam
        samtools index ${X}.bam
        """
}

process SamtoolsMerge{
    //  samtools merge – merges multiple sorted files into a single file 
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'biocontainers/samtools:v1.7.0_cv4'
    label 'many_cpu_medium'

    input:
        path(input_sam)
        val output_base

    output:
        path "${output_base}.bam"

    script:
        """
        samtools merge -O bam ${output_base}.bam $input_sam
        """
}

process SamtoolsQuickcheck{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'biocontainers/samtools:v1.7.0_cv4'
    label 'many_cpu_medium'

    input:
        path(input_sam)

    output:
        stdout

    script:
        """
        samtools quickcheck $input_sam && echo 'Samtools quickcheck ok' || echo 'Samtools quickcheck fail!'
        """
}

process SamtoolsFlagstats{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'biocontainers/samtools:v1.7.0_cv4'
    label 'many_cpu_medium'

    input:
        path(input_sam)
    
    output:
        stdout

    script:
        """
        samtools flagstat $input_sam
        """
}

process SamtoolsFlagstatsMapPercentage{
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'biocontainers/samtools:v1.7.0_cv4'
    label 'many_cpu_medium'

    input:
        tuple val(X), path(bam), path(bai)
    
    output:
        tuple val(X), stdout

    script:
    //         samtools flagstat $bam | grep % | cut -d '(' -f 2 | cut -d '%' -f1

        """
        set -e -o pipefail
        samtools flagstat $bam | grep % | cut -d '(' -f 2 | cut -d '%' -f1 | head -n 1
        """ 
}

process SamtoolsMergeTuple{
    //  merge n number of bams into one
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'biocontainers/samtools:v1.7.0_cv4'
    label 'many_cpu_medium'

    input:
        tuple val(X), path(bam), path(bai)

    output:
        tuple val(X), path("${X}.merged.bam"), path("${X}.merged.bam.bai")
    
    script:
    // contains bug, remove in future release
    println('Warning Depreciated')
    """
    samtools merge -O bam ${X}.merged.bam \$(find . -name '*.bam')
    samtools index ${X}.merged.bam
    """
}

process SamtoolsMergeBams{
    //  merge n number of bams into one
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'biocontainers/samtools:v1.7.0_cv4'
    label 'many_cpu_medium'

    input:
        val(X)
        file(bam_in)

    output:
        tuple val(X), path("${X}.merged.bam"), path("${X}.merged.bam.bai")
    
    script:
    """
    ls
    samtools merge -O bam ${X}.merged.bam \$(find . -name '*.bam')
    samtools index ${X}.merged.bam
    """
}

process RemoveUnmappedReads{
    // Sort, convert and index 
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'biocontainers/samtools:v1.7.0_cv4'
    label 'many_cpu_medium'

    input:
        tuple val(X), path(bam_in), path(bai_in)

    output:
        tuple val(X), path("${X}.mapped.bam"), path("${X}.mapped.bam.bai")

    script:
        """
        samtools view -b -F 4 $bam_in > ${X}.mapped.bam
        samtools index ${X}.mapped.bam
        """
}

process PrimaryMappedFilter{
    // Sort, convert and index 
    publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'biocontainers/samtools:v1.7.0_cv4'
    label 'many_cpu_medium'

    input:
        tuple val(X), path(bam_in), path(bai_in)

    output:
        tuple val(X), path("${X}.primary_mapped.bam"), path("${X}.primary_mapped.bam.bai")

    script:
        """
        samtools view -b -F 256 $bam_in > ${X}.primary_mapped.bam
        samtools index ${X}.primary_mapped.bam
        """
}
