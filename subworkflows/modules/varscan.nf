#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process VarscanFiltered {
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        tuple val(X), path(input_bam_file),path(input_bai_file), path(reference)

    output:
        tuple val(X), path("${input_bam_file.SimpleName}.vcf")
        // tuple val(X), path("${input_bam_file.SimpleName}.vcf"), path("${X}_both.vcf"), path("${X}_fwd_rev_support_snps.vcf"), path("${X}_rev_snps.vcf"), path("${X}_fwd_snps.vcf")

    script:
    // TODO: filter on depth and make empty vcf if depth not reached with warning to stdout
    """
        # -r can be used to limit the variant calling in mpileup
        # Call using only reverse reads
        samtools mpileup --incl-flags 0x10 -q ${params.varscan.min_mq} -Q ${params.varscan.min_bqual} -B -d ${params.varscan.max_depth} -A -f $reference $input_bam_file | varscan mpileup2cns --variants 1 --output-vcf 1 --min-coverage ${params.varscan.min_coverge} --min-avg-qual ${params.varscan.min_avg_qual} --min-var-freq ${params.varscan.min_var_freq} --strand-filter 0 --min-reads2 ${params.varscan.min_support_reads} | bcftools view -v snps > ${X}_rev_snps.vcf
        # Call using only fwd reads
        samtools mpileup --excl-flags 0x714 -q ${params.varscan.min_mq} -Q ${params.varscan.min_bqual} -B -d ${params.varscan.max_depth} -A -f $reference $input_bam_file | varscan mpileup2cns --variants 1 --output-vcf 1 --min-coverage ${params.varscan.min_coverge} --min-avg-qual ${params.varscan.min_avg_qual} --min-var-freq ${params.varscan.min_var_freq} --strand-filter 0 --min-reads2 ${params.varscan.min_support_reads} | bcftools view -v snps > ${X}_fwd_snps.vcf
        # Find calls that are in both
        bedtools intersect -header -u -a ${X}_rev_snps.vcf -b ${X}_fwd_snps.vcf > ${X}_fwd_rev_support_snps.vcf

        # Call indels using only reverse reads
        samtools mpileup --incl-flags 0x10 -q ${params.varscan.min_mq} -Q ${params.varscan.min_bqual} -B -d ${params.varscan.max_depth} -A -f $reference $input_bam_file | varscan mpileup2cns --variants 1 --output-vcf 1 --min-coverage ${params.varscan.min_coverge} --min-avg-qual ${params.varscan.min_avg_qual} --min-var-freq ${params.varscan.min_indel_freq} --strand-filter 0 --min-reads2 ${params.varscan.min_support_reads} | bcftools view -v indels > ${X}_rev_indels.vcf
        # Call using only fwd reads
        samtools mpileup --excl-flags 0x714 -q ${params.varscan.min_mq} -Q ${params.varscan.min_bqual} -B -d ${params.varscan.max_depth} -A -f $reference $input_bam_file | varscan mpileup2cns --variants 1 --output-vcf 1 --min-coverage ${params.varscan.min_coverge} --min-avg-qual ${params.varscan.min_avg_qual} --min-var-freq ${params.varscan.min_indel_freq} --strand-filter 0 --min-reads2 ${params.varscan.min_support_reads} | bcftools view -v indels > ${X}_fwd_indels.vcf
        # Find calls that are in both
        bedtools intersect -header -u -a ${X}_rev_indels.vcf -b ${X}_fwd_indels.vcf > ${X}_fwd_rev_support_indels.vcf

        # Call using forward and reverse
        samtools mpileup -q ${params.varscan.min_mq} -Q ${params.varscan.min_bqual} -B -d ${params.varscan.max_depth} -A -f $reference $input_bam_file | varscan mpileup2cns --variants 1 --output-vcf 1 --min-coverage ${params.varscan.min_coverge} --min-avg-qual ${params.varscan.min_avg_qual} --min-var-freq ${params.varscan.min_var_freq} --strand-filter 1 --min-reads2 ${params.varscan.min_support_reads} > ${X}_both.vcf

        # Keep only calls that were supported by both strands and remove indels that have a smaller AF than min_indel_freq
        bedtools intersect -header -u -a ${X}_both.vcf -b ${X}_fwd_rev_support_snps.vcf  > ${X}_filtered_snps.vcf
        bedtools intersect -header -u -a ${X}_both.vcf -b ${X}_fwd_rev_support_indels.vcf  > ${X}_filtered_indels.vcf

        bcftools view ${X}_filtered_snps.vcf > ${input_bam_file.SimpleName}.vcf
    """
}
