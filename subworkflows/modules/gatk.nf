#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process CreateRererenceDict {
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'broadinstitute/gatk:4.2.3.0'
    label 'many_med_cpu_huge_mem'


    input:
        path(reference)

    output:
        path("${reference.simpleName}.dict")

    script:
        """
        gatk CreateSequenceDictionary -R $reference
        """
} 

process GatkAddOrReplaceReadGroups {
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'broadinstitute/gatk:4.2.3.0'
    label 'many_med_cpu_huge_mem'


    input:
        tuple val(X), path(input_bam_file), path(input_bai_file)

    output:
        tuple val(X), path("${X}.RG.bam"), path("${X}.RG.bam.bai")


    script:
        """
        gatk AddOrReplaceReadGroups \
        I=$input_bam_file \
        O=${X}.RG.bam \
        RGID=4 \
        RGLB=lib1 \
        RGPL=ILLUMINA \
        RGPU=unit1 \
        RGSM=20
        samtools index ${X}.RG.bam
        """
}

process Mutect2TumorOnly {
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'broadinstitute/gatk:4.2.3.0'
    label 'few_cpu_intensive'

    input:
        tuple val(X), path(input_bam), path(input_bai), path(reference), path(reference_dict)

    output:
        tuple val(X), path("${X}.mutect2.vcf")

    script:
        """
        samtools faidx $reference
        gatk Mutect2 \
        --native-pair-hmm-threads ${task.cpus} \
        -R $reference \
        -I $input_bam \
        -O ${X}.mutect2.vcf
    """
}