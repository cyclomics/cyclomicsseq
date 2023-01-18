


process AddDepthToJson{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    label 'many_cpu_medium'

    input:
        tuple val(X), path(tidehuntertable)
        path(depth_json)

    output:
        tuple val(X), path("${X}.depth.json")

    script:
        """  
        add_depth_info_json.py --global_json $tidehuntertable --depth_json $depth_json --output ${X}.depth.json
        """
}

process AnnotateBamXTags{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/consensus_aligned", mode: 'copy'

    label 'many_low_cpu_high_mem'

    input:
        tuple val(X), path(bam), path(bai)
        path(sequencing_summary)

    output:
        tuple val(X), path("${X}.taged.bam"), path("${X}.taged.bam.bai")

    script:
        """
        annotate_bam_x.py $sequencing_summary $bam ${X}.taged.bam
        """
}

process AnnotateBamYTags{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'

    label 'many_low_cpu_high_mem'

    input:
        tuple val(X), path(bam), path(json)

    output:
        tuple val(X), path("${X}.annotated.bam"), path("${X}.annotated.bam.bai")

    script:
        """
        samtools index $bam
        annotate_bam_y.py $json $bam ${X}.annotated.bam
        """
}

process CollectClassificationTypes{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'

    input:
        path(metadata_json)

    output:
        path("classification_count.txt")
    
    script:
        """
        gather_readtypes.py "*.metadata.json" classification_count.txt
        """
}

process FindVariants{
    publishDir "${params.output_dir}/variants", mode: 'copy'

    label 'max_performance'

    input:
        path(reference_genome)
        tuple val(X), path(bam), path(bai)
        tuple val(X), path(validation_bed)

    output:
        tuple path("${bam.simpleName}.snp.vcf"), path("${bam.simpleName}.indel.vcf")
    
    script:
        """
        determine_vaf.py $reference_genome $validation_bed $bam ${bam.simpleName}.snp.vcf ${bam.simpleName}.indel.vcf --threads ${task.cpus} 
        """
}

process FilterVariants{
    publishDir "${params.output_dir}/variants", mode: 'copy'

    input:
        tuple path(snp_vcf), path(indel_vcf), val(X), path(perbase_table)

    output:
        tuple path("${snp_vcf.simpleName}_filtered.snp.vcf"), path("${snp_vcf.simpleName}_filtered.indel.vcf")
    
    script:
        """
        vcf_filter.py -i $snp_vcf -o ${snp_vcf.simpleName}_filtered.snp.vcf -p $perbase_table \
        --min_dir_ratio $params.snp_filters.min_dir_ratio \
        --min_dir_count $params.snp_filters.min_dir_count \
        --min_dpq $params.snp_filters.min_dpq \
        --min_dpq_n $params.snp_filters.min_dpq_n \
        --min_dpq_ratio $params.snp_filters.min_dpq_ratio \
        --min_vaf $params.snp_filters.min_vaf \
        --min_rel_ratio $params.snp_filters.min_rel_ratio \
        --min_abq $params.snp_filters.min_abq

        vcf_filter.py -i $indel_vcf -o ${indel_vcf.simpleName}_filtered.indel.vcf -p $perbase_table \
        --min_dir_ratio $params.indel_filters.min_dir_ratio \
        --min_dir_count $params.indel_filters.min_dir_count \
        --min_dpq $params.indel_filters.min_dpq \
        --min_dpq_n $params.indel_filters.min_dpq_n \
        --min_dpq_ratio $params.indel_filters.min_dpq_ratio \
        --min_vaf $params.indel_filters.min_vaf \
        --min_rel_ratio $params.indel_filters.min_rel_ratio \
        --min_abq $params.indel_filters.min_abq
        """
}

process MergeNoisyVCF{
    publishDir "${params.output_dir}/variants", mode: 'copy'

    input:
        tuple path(noisy_snp_vcf), path(noisy_indel_vcf)

    output:
        path("${noisy_snp_vcf.simpleName}.noisy_merged.vcf")
    
    script:
        """
        bcftools sort $noisy_snp_vcf -o ${noisy_snp_vcf.simpleName}.sorted.snp.vcf
        bgzip ${noisy_snp_vcf.simpleName}.sorted.snp.vcf
        tabix ${noisy_snp_vcf.simpleName}.sorted.snp.vcf.gz

        bcftools sort $noisy_indel_vcf -o ${noisy_indel_vcf.simpleName}.sorted.indel.vcf
        bgzip ${noisy_indel_vcf.simpleName}.sorted.indel.vcf
        tabix ${noisy_indel_vcf.simpleName}.sorted.indel.vcf.gz

        bcftools concat -a ${noisy_snp_vcf.simpleName}.sorted.snp.vcf.gz ${noisy_indel_vcf.simpleName}.sorted.indel.vcf.gz \
        -O v -o ${noisy_snp_vcf.simpleName}.noisy_merged.tmp.vcf
        bcftools sort ${noisy_snp_vcf.simpleName}.noisy_merged.tmp.vcf -o ${noisy_snp_vcf.simpleName}.noisy_merged.vcf
        rm ${noisy_snp_vcf.simpleName}.noisy_merged.tmp.vcf
        """
}

process MergeFilteredVCF{
    publishDir "${params.output_dir}/variants", mode: 'copy'

    input:
        tuple path(filtered_snp_vcf), path(filtered_indel_vcf)

    output:
        path("${filtered_snp_vcf.simpleName}.filtered_merged.vcf")
    
    script:
        """
        bcftools sort $filtered_snp_vcf -o ${filtered_snp_vcf.simpleName}.sorted.snp.vcf
        bgzip ${filtered_snp_vcf.simpleName}.sorted.snp.vcf
        tabix ${filtered_snp_vcf.simpleName}.sorted.snp.vcf.gz

        bcftools sort $filtered_indel_vcf -o ${filtered_indel_vcf.simpleName}.sorted.indel.vcf
        bgzip ${filtered_indel_vcf.simpleName}.sorted.indel.vcf
        tabix ${filtered_indel_vcf.simpleName}.sorted.indel.vcf.gz

        bcftools concat -a ${filtered_snp_vcf.simpleName}.sorted.snp.vcf.gz ${filtered_indel_vcf.simpleName}.sorted.indel.vcf.gz \
        -O v -o ${filtered_snp_vcf.simpleName}.filtered_merged.tmp.vcf
        bcftools sort ${filtered_snp_vcf.simpleName}.filtered_merged.tmp.vcf -o ${filtered_snp_vcf.simpleName}.filtered_merged.vcf
        rm ${filtered_snp_vcf.simpleName}.filtered_merged.tmp.vcf
        """
}

process AnnotateVCF{
    publishDir "${params.output_dir}/variants", mode: 'copy'

    input:
        path(variant_vcf)


    output:
        path("${variant_vcf.simpleName}_annotated.vcf")
    
    script:
        """
        annotate_vcf.py $variant_vcf ${variant_vcf.simpleName}_annotated.vcf
        """
}

process PlotFastqsQUalAndLength{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/QC", mode: 'copy'

    input:
        path(fastq)
        val grep_pattern
        val plot_file_prefix
        val tab_name

    output:
        tuple path("${plot_file_prefix}_histograms.html"), path("${plot_file_prefix}_histograms.json")
    
    script:
        """
        plot_fastq_histograms.py $grep_pattern ${plot_file_prefix}_histograms.html $tab_name
        """
}

process PlotReadStructure{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/QC", mode: 'copy'

    input:
        tuple val(X), path(bam), path(bai)

    output:
        tuple path("${bam.simpleName}_aligned_segments.html"), path("${bam.simpleName}_read_structure.html"), path("${bam.simpleName}_read_structure.json")

    script:
        """
        plot_read_structure_donut.py $bam ${bam.simpleName}_aligned_segments.html ${bam.simpleName}_read_structure.html
        """
}

process PlotVcf{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/QC", mode: 'copy'

    input:
        path(vcf)

    output:
        tuple path("${vcf.simpleName}.html"), path("${vcf.simpleName}.json")
    
    script:
        """
        plot_vcf.py $vcf ${vcf.simpleName}.html

        """
}

process PasteVariantTable{
    publishDir "${params.output_dir}/QC", mode: 'copy'

    input:
        path(vcf_file)

    output:
        path("${vcf_file.simpleName}_table.json")
    
    script:
        """
        write_variants_table.py $vcf_file ${vcf_file.simpleName}_table.json 'Variant Table'

        """
}

process PlotQScores{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/QC", mode: 'copy'

    input:
        tuple val(X), path(split_pileup)
        tuple val(Y), path(consensus_pileup)

    output:
        tuple path("${consensus_pileup.simpleName}.html"), path("${consensus_pileup.simpleName}.csv"), path("${consensus_pileup.simpleName}.json")
    
    script:
        """
        plot_bam_accuracy.py $split_pileup $consensus_pileup ${consensus_pileup.simpleName}.html
        """
}

process PlotMetadataStats{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/QC", mode: 'copy'

    input:
        path(jsons)

    output:
        tuple path("metadata_plots.html"), path("metadata_plots.json")
    
    script:
        """
        plot_metadata.py . metadata_plots.html
        """
}
process PlotReport{
    publishDir "${params.output_dir}/QC", mode: 'copy'

    input:
        path(jsons)

    output:
        path("report.html")
    
    script:
        """
        generate_report.py '${params}' $workflow.manifest.version
        """

}
