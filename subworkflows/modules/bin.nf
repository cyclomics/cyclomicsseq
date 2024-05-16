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
    // publishDir "${params.output_dir}/consensus_aligned_tagged", mode: 'copy'
    label 'few_very_memory_intensive'

    input:
        tuple val(X), path(bam), path(bai)
        path(sequencing_summary)

    output:
        tuple val(X), path("${X}.tagged.bam"), path("${X}.tagged.bam.bai")

    script:
        """
        annotate_bam_x.py $sequencing_summary $bam ${X}.tagged.bam
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
    label 'many_low_cpu_high_mem'

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
    // publishDir "${params.output_dir}/variants", mode: 'copy'
    label 'max_performance'

    input:
        path(reference_genome)
        tuple val(X), path(bam), path(bai)
        tuple val(X), path(validation_bed)

    output:
        tuple val(X), path("${bam.simpleName}_snp.vcf"), path("${bam.simpleName}_indel.vcf")
    
    script:
        // We sleep and access the reference genome, since in some rare cases the file needs accessing to 
        // not cause issues in the python code. 
        """
        sleep 1
        ls
        head $reference_genome
        determine_vaf.py $reference_genome $validation_bed $bam ${bam.simpleName}_snp.vcf ${bam.simpleName}_indel.vcf --threads ${task.cpus} 
        """
}

process FilterValidateVariants{
    // publishDir "${params.output_dir}/variants", mode: 'copy'
    label 'many_low_cpu_high_mem'

    input:
        tuple val(X), path(vcf), path(perbase_table)
        val(mode)

    output:
        tuple val(X), path("${vcf.simpleName}_filtered.vcf")
    
    script:
    if ( mode == 'snp' )
        """
        vcf_filter.py -i $vcf -o ${vcf.simpleName}_filtered.vcf \
        --perbase_table $perbase_table \
        --dynamic_vaf_params $params.dynamic_vaf_params_file \
        --min_ao $params.snp_filters.min_ao \
        --min_dpq $params.snp_filters.min_dpq \
        --min_dpq_n $params.snp_filters.min_dpq_n \
        --min_dpq_ratio $params.snp_filters.min_dpq_ratio \
        --max_sap 0 \
        --min_rel_ratio $params.snp_filters.min_rel_ratio \
        --min_abq $params.snp_filters.min_abq
        """

    else if ( mode == 'indel' )
        """
        vcf_filter.py -i $vcf -o ${vcf.simpleName}_filtered.vcf \
        --perbase_table $perbase_table \
        --dynamic_vaf_params $params.dynamic_vaf_params_file \
        --min_ao $params.indel_filters.min_ao \
        --min_dpq $params.indel_filters.min_dpq \
        --min_dpq_n $params.indel_filters.min_dpq_n \
        --min_dpq_ratio $params.indel_filters.min_dpq_ratio \
        --max_sap 0 \
        --min_rel_ratio $params.indel_filters.min_rel_ratio \
        --min_abq $params.indel_filters.min_abq   
        """

    else
        error "Invalid merging mode (VCF type): ${mode}. Valid modes are 'snp' or 'indel'."
}

process SortVCF{
    // publishDir "${params.output_dir}/variants", mode: 'copy'
    label 'many_low_cpu_high_mem'

    input:
        tuple val(X), path(vcf)

    output:
        tuple val(X), path("${vcf.simpleName}_sorted.vcf.gz"), path("${vcf.simpleName}_sorted.vcf.gz.tbi")
    
    script:
        """
        bcftools sort $vcf -o ${vcf.simpleName}_sorted.vcf --temp-dir .
        bgzip -f ${vcf.simpleName}_sorted.vcf
        tabix -f ${vcf.simpleName}_sorted.vcf.gz
        """
}


process MergeVCF{
    publishDir "${params.output_dir}/variants", mode: 'copy'
    label 'many_low_cpu_high_mem'

    input:
        tuple val(X), path(snp_vcf), path(snp_tbi), path(indel_vcf), path(indel_tbi)
        val(mode)

    output:
        tuple val(X), path("${X}_${mode}.vcf")
    
    script:
    if ( mode == 'noisy' )
        """
        bcftools concat -a ${snp_vcf.simpleName}.vcf.gz ${indel_vcf.simpleName}.vcf.gz \
        -O v -o ${X}_${mode}.tmp.vcf \
        2> >(tee -a error.txt >&2) && \
        if grep -q "E::" error.txt; then \
        echo "VCF concatenation failed in certain rows.\n"; \
        exit 1; \
        fi

        bcftools sort ${X}_${mode}.tmp.vcf -o ${X}_${mode}.vcf --temp-dir .
        rm ${X}_${mode}.tmp.vcf
        """

    else if ( mode == 'filtered' )
        """
        bcftools concat -a ${snp_vcf.simpleName}.vcf.gz ${indel_vcf.simpleName}.vcf.gz \
        -O v -o ${X}_${mode}.tmp.vcf \
        2> >(tee -a error.txt >&2) && \
        if grep -q "E::" error.txt; then \
        echo "VCF concatenation failed in certain rows.\n"; \
        exit 1; \
        fi

        bcftools sort ${X}_${mode}.tmp.vcf -o ${X}_${mode}.vcf --temp-dir .
        rm ${X}_${mode}.tmp.vcf
        """

    else
        error "Invalid merging mode (VCF type): ${mode}. Valid modes are 'noisy' or 'filtered'."
}


process AnnotateVCF{
    publishDir "${params.output_dir}/variants", mode: 'copy'
    label 'many_low_cpu_high_mem'

    input:
        tuple val(X), path(vcf)


    output:
        tuple val(X), path("${vcf.simpleName}_annotated.vcf")
    
    script:
        """
        annotate_vcf.py $vcf ${vcf.simpleName}_annotated.vcf
        """
}

process PlotFastqsQUalAndLength{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/QC", mode: 'copy'
    label 'many_low_cpu_high_mem'
    errorStrategy 'ignore'

    input:
        path(fastq)
        val grep_pattern
        val plot_file_prefix
        val tab_name

    output:
        path("${plot_file_prefix}_histograms.json")
    
    script:
        """
        plot_fastq_histograms.py $grep_pattern ${plot_file_prefix}_histograms.html $tab_name $params.priority_limit \
        2> >(tee -a error.txt >&2) || catch_plotting_errors.sh error.txt ${plot_file_prefix}_histograms.json
        """
}

process PlotReadStructure{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/QC", mode: 'copy'
    label 'many_cpu_medium'

    memory { task.memory * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 3

    input:
        tuple val(X), path(bam), path(bai)

    output:
        path("${bam.simpleName}_read_structure.json")

    script:
        // This takes a lot of RAM when the sequencing summary is big!
        """
        samtools sort -n -o tmp_readname_sorted_${bam.simpleName}.bam ${bam}
        plot_read_structure_donut.py tmp_readname_sorted_${bam.simpleName}.bam \
        ${bam.simpleName}_read_structure.json \
        $params.priority_limit \
        2> >(tee -a error.txt >&2) || catch_plotting_errors.sh error.txt ${bam.simpleName}_read_structure.json
        """
}

process PlotVcf{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/QC", mode: 'copy'
    label 'many_low_cpu_high_mem'
    memory { task.memory * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 3

    input:
        path(vcf)

    output:
        path("${vcf.simpleName}.json")
    
    script:
        // This takes a lot of RAM when the sequencing summary is big!
        """
        plot_vcf.py $vcf ${vcf.simpleName}.html $params.priority_limit \
        2> >(tee -a error.txt >&2) || catch_plotting_errors.sh error.txt ${vcf.simpleName}.json
        """
}

process PasteVariantTable{
    publishDir "${params.output_dir}/QC", mode: 'copy'
    label 'many_low_cpu_high_mem'
    errorStrategy 'ignore'

    input:
        path(vcf_file)

    output:
        path("${vcf_file.simpleName}_table.json")
    
    script:
        """
        write_variants_table.py $vcf_file ${vcf_file.simpleName}_table.json 'Variant table' $params.priority_limit \
        2> >(tee -a error.txt >&2) || catch_plotting_errors.sh error.txt ${vcf_file.simpleName}_table.json
        """
}

process PlotQScores{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/QC", mode: 'copy'
    label 'many_low_cpu_high_mem'
    errorStrategy 'ignore'

    input:
        tuple val(X), path(split_pileup)
        tuple val(Y), path(consensus_pileup)

    output:
        path("${consensus_pileup.simpleName}.json")
    
    script:
        """
        plot_bam_accuracy.py $split_pileup $consensus_pileup ${consensus_pileup.simpleName}.html $params.priority_limit \
        2> >(tee -a error.txt >&2) || catch_plotting_errors.sh error.txt ${consensus_pileup.simpleName}.json
        """
}

process PlotMetadataStats{
    // publishDir "${params.output_dir}/${task.process.replaceAll(':', '/')}", pattern: "", mode: 'copy'
    publishDir "${params.output_dir}/QC", mode: 'copy'
    label 'many_low_cpu_huge_mem'
    memory { task.memory * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 3

    input:
        path(jsons)

    output:
        path("metadata_plots.json")
    
    script:
        """
        plot_metadata.py . metadata_plots.html $params.priority_limit --subsample_size $params.metadata.subsample_size \
        2> >(tee -a error.txt >&2) || catch_plotting_errors.sh error.txt metadata_plots.json
        """
}

process PlotReport{
    publishDir "${params.output_dir}", mode: 'copy'
    label 'many_low_cpu_high_mem'

    input:
        path(jsons)

    output:
        path("report.html")
    
    script:
        """
        generate_report.py '${params}' $workflow.manifest.version $params.priority_limit
        """
}