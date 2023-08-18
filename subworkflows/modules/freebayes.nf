#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process Freebayes {
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    container 'nfcore/sarek:2.7.1'
    label 'many_cpu_medium'

    input:
        tuple val(X), path(input_bam_file), path(input_bai_file), path(reference)
        tuple val(X), path(roi)

    output:
        path("${input_bam_file.SimpleName}.vcf")

    script:
        ref = reference
        """
        freebayes \
        --min-alternate-fraction $params.min_alt_fraction
        --min-alternate-count $params.min_alt_count
        --min-base-quality $params.min_base_qual
        -t $roi \
        -f $ref \
        --vcf ${input_bam_file.SimpleName}.vcf \
        $input_bam_file
        """
}

process FilterFreebayesVariants{
    publishDir "${params.output_dir}/variants", mode: 'copy'
    label 'many_low_cpu_high_mem'

    input:
        tuple path(vcf), val(X), path(perbase_table)

    output:
        tuple path("${snp_vcf.simpleName}_filtered.snp.vcf"), path("${snp_vcf.simpleName}_filtered.indel.vcf")
    
    script:
        """
        vcf_filter.py -i $vcf -o ${snp_vcf.simpleName}_filtered.snp.vcf \
        --perbase_table $perbase_table \
        --dynamic_vaf_params $params.dynamic_vaf_params_file \
        --min_dpq $params.var_filters.min_dpq \
        --min_dpq_n $params.var_filters.min_dpq_n \
        --min_dpq_ratio $params.var_filters.min_dpq_ratio \
        --max_sap $params.var_filters.max_sap \
        --min_rel_ratio $params.var_filters.min_rel_ratio \
        --min_abq $params.var_filters.min_abq        
        """
}